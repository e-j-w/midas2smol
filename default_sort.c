#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "config.h"
#include "grif-format.h"
#include "smol-format.h"
#include "default_sort.h"

int          odb_daqsize;// number of daq channels currently defined in the odb
int         subsys_table[MAX_DAQSIZE];
int        crystal_table[MAX_DAQSIZE]; // Ge/BGO have 4 "crystals" per clover
int        element_table[MAX_DAQSIZE]; // BGO have 5 elements per crystal
int       polarity_table[MAX_DAQSIZE]; // 1 is negative, 0 is positive, -1 is unset
short       address_chan[MAX_ADDRESS];
static short  addr_table[MAX_DAQSIZE]; short   *addrs = addr_table;
       char    chan_name[MAX_DAQSIZE][CHAN_NAMELEN];
static int   dtype_table[MAX_DAQSIZE]; int    *dtypes = dtype_table;
static float  gain_table[MAX_DAQSIZE]; float   *gains = gain_table;
static float  offs_table[MAX_DAQSIZE]; float *offsets = offs_table;
static float  quad_table[MAX_DAQSIZE]; float   *quads = quad_table;
float  pileupk1[MAX_DAQSIZE][7];
float  pileupk2[MAX_DAQSIZE][7];
float  pileupE1[MAX_DAQSIZE][7];
float  crosstalk[MAX_DAQSIZE][3][16];
static short *chan_address = addr_table;
static int subsys_initialized[MAX_SUBSYS];
extern Grif_event grif_event[PTR_BUFSIZE];
extern uint64_t psd_vals[MAX_PSD_VALS];

int presort_window_width= 1940;  // 19.4us needed for all crosstalk corrections. 5us needed for pileup corrections.
int sort_window_width   = 200; //  2us - MAXIMUM (indiv. gates can be smaller)


// Default sort function declarations
extern int perform_pileup_correction(Grif_event *ptr, Grif_event *alt, int dt, int chan, int chan2, int i, int end_idx);
extern uint8_t fill_smol_entry(FILE *out, const int win_idx, const int frag_idx);

//generates a random double value on the interval [0,1]
//for smoothing purposes
double randomDbl(){
  return ((double)(rand()) / (double)(RAND_MAX));
}

//get the CFD corrected time
//implemented in GRSISort at https://github.com/GRIFFINCollaboration/GRSIData/blob/4b5dbe964d18190c03e151bca818d980c03a9bfc/libraries/TGRSIFormat/TGRSIMnemonic.cxx#L218
double getGrifTime(Grif_event *ptr){
  /*int cfdTSDiff = abs((int)(ptr->ts & 0x3ffff) - (int)(ptr->cfd >> 4));
  if(cfdTSDiff >= 30){
    return -1.0; //discard hit, failed CFD
  }*/
  double timeNs = ((unsigned long)(ptr->ts & ~0x000000000003ffffUL))*10.0; //timestamp value, excluding lower 18 bits
  //printf("timeNs: %f\n",timeNs);
  /*if(timeNs > 1.0E19){
    printf("getGrifTime - ts: %lu, time in ns: %f",ptr->ts, timeNs);
    getc(stdin);
  }*/
  return (timeNs + (((double)ptr->cfd + randomDbl())/1.6));
  //return ((double)ptr->ts + randomDbl())*10.0;
}

// odb tables need to be transferred into config, which is saved with histos
int init_default_sort(Config *cfg, Sort_status *arg)
{
  Cal_coeff *cal;
  int i, j;

  // Initialize all pileup parameters to unset values
  for(i=0; i<odb_daqsize; i++){
    for(j=0; j<7; j++){
      pileupk1[i][j] = pileupk2[i][j] = pileupE1[i][j] = -1;
    }
    for(j=0; j<16; j++){
      crosstalk[i][0][j] = crosstalk[i][1][j] = crosstalk[i][2][j] = -1;
    }
  }

  cfg->odb_daqsize = odb_daqsize;
  for(i=0; i<odb_daqsize; i++){ // update config with odb info
    edit_calibration(cfg, chan_name[i], offsets[i], gains[i], quads[i], pileupk1[i], pileupk2[i], pileupE1[i],
      crosstalk[i][0], crosstalk[i][1], crosstalk[i][2], chan_address[i],  dtype_table[i], arg->cal_overwrite );
  }
  // ALSO need to transfer config info to the arrays that are used in sort
  for(i=0; i<odb_daqsize; i++){

    cal = cfg->calib[i];
    if( strcmp(chan_name[i], cal->name) != 0 ){ // conf not in odb order
      for(j=0; j<cfg->ncal; j++){ cal = cfg->calib[j];
        if( strcmp(chan_name[i], cal->name) == 0 ){ break; }
      }
      if( j == cfg->ncal ){ continue; } // not found in config
    }

    // overwrite = 0 => USE CONFIG NOT ODB for offset, gain, quads
    if( arg->cal_overwrite == 0 ){
      offsets[i]=cal->offset; gains[i]=cal->gain;  quads[i]=cal->quad;
    }

    // Pileup parameters do not exist in the MIDAS ODB so must always be copied from the config
    for(j=0; j<7; j++){
      pileupk1[i][j] = (isnan(cal->pileupk1[j])) ? 0.0 : cal->pileupk1[j];
      pileupk2[i][j] = (isnan(cal->pileupk2[j])) ? 0.0 : cal->pileupk2[j];
      pileupE1[i][j] = (isnan(cal->pileupE1[j])) ? 0.0 : cal->pileupE1[j];
    }
    for(j=0; j<16; j++){
      crosstalk[i][0][j] = (isnan(cal->crosstalk0[j])) ? 0.0 : cal->crosstalk0[j];
      crosstalk[i][1][j] = (isnan(cal->crosstalk1[j])) ? 0.0 : cal->crosstalk1[j];
      crosstalk[i][2][j] = (isnan(cal->crosstalk2[j])) ? 0.0 : cal->crosstalk2[j];
    }
  }

  return(0);
}

//#######################################################################
//#####        BASIC DEFAULT SORT (common to most experiments)      #####
//#######################################################################

float spread(int val){ return( val + rand()/(1.0*RAND_MAX) ); }
int GetIDfromAddress(unsigned short addr){ // address must be an unsigned short
  return(address_chan[addr]);
}


//#######################################################################
//######## PRESORT(gain corrections, suppression)              ##########
//#######################################################################

// We don't do addback here, instead passing through all single crystal hits
// to the output file. Addback can be implemented when sorting the output
// file(s).

// (used to be called apply_gains) this is the first function to be called
// on processing an event - before any singles/coinc-sorting ...
// ** the current event has just been added and is last in the window
//      => all other window events are BEFORE the current event
int pre_sort_enter(int start_idx, int frag_idx)
{
  Grif_event *alt, *ptr = &grif_event[frag_idx];
  float energy, correction;
  int dt, bin, chan2, chan = ptr->chan;
  int clover, ge1, c1,c2, add;
  int ct_index[4][4] = {{-1,0,1,2},{0,-1,1,2},{0,1,-1,2},{0,1,2,-1}};

  // Protect against invalid channel numbers
  if( chan < 0 || chan >= odb_daqsize ){
    if( ptr->address == 0xFFFF ){
      /*
      ppg_index=-1;
      for(i=0; i<N_PPG_PATTERNS; i++){ if( (ptr->master_pattern & 0xFFFF) == ppg_patterns[i] ){ ppg_index = i; break; } }
      if(ppg_index<0){ fprintf(stderr,"unrecognized ppg pattern, 0x%04X\n", (ptr->master_pattern & 0xFFFF)); return(-1); }
      */
      //  fprintf(stdout,"PPG PATTERN: 0x%04X (%s, %s) at time %10.4f seconds\n", (ptr->master_pattern & 0xFFFF), ppg_handles[ppg_index], ppg_names[ppg_index], (double)(ptr->ts/100000000) );
    } else {
      fprintf(stderr,"unpack_event: ignored event in chan:%d [0x%04x]\n", chan, ptr->address );
    }
    return(-1);
  }

  // Calculate the energy and calibrated energies
  energy = ( ptr->integ1 == 0 ) ? ptr->q1 : spread(ptr->q1)/ptr->integ1;
  ptr->ecal = ptr->esum=offsets[chan]+energy*(gains[chan]+energy*quads[chan]);
  // NOBODY CURRENTLY USES e2,e3,e4 ...

  // Assign the subsys type
  if( (ptr->subsys = subsys_table[chan]) == -1 ){ return(-1); }

    // HPGe B
    if( ptr->subsys == SUBSYS_HPGE_B){
      ptr->pu_class = PU_OTHER; // Pileup class - default value for all HPGe events
      if(ptr->pileup==1 && ptr->nhit ==1){
        ptr->pu_class = PU_SINGLE_HIT; // Single hit events, no pileup, this is the most common type of HPGe event
      }
    }

    // HPGe A
    if( ptr->subsys == SUBSYS_HPGE_A){
      ptr->pu_class = PU_OTHER; // Pileup class - default value for all HPGe events
      if(ptr->pileup==1 && ptr->nhit ==1){
        ptr->pu_class = PU_SINGLE_HIT; // Single hit events, no pileup, this is the most common type of HPGe event
      }

      // HPGe Clover time-dependant crosstalk corrections within same clover
      int i = start_idx;
      while( i != frag_idx ){ // need at least two events in window
        if( ++i >=  PTR_BUFSIZE ){ i=0; } alt = &grif_event[i]; // WRAP
        if( (chan2 = alt->chan)<0 || alt->chan >= odb_daqsize ){
          fprintf(stderr,"presort error: ignored event in chan:%d\n",alt->chan );
          continue;
        }
        if((dt=ptr->ts - alt->ts)>479 || alt->subsys != SUBSYS_HPGE_A){ continue; }

        if(chan2 != chan ){
          if((clover=(int)(crystal_table[chan2]/4)) == (int)(crystal_table[chan]/4)){
            // HPGe Clover time-dependant crosstalk corrections within same clover
            // dt is always positive here
            // The original hit (ptr) came after the crosstalk-inducing hit (alt)
            // Make correction to ptr hit based on energy of alt.
            bin = (int)((1940+dt)/160);
            if(bin<0 || bin>15){ fprintf(stderr,"pre_sort_enter bin [%d] out of bounds for dt %d\n",bin,dt); continue; }
            ge1 = crystal_table[chan];
            c1 = ge1%4;
            c2 = ct_index[c1][(int)(crystal_table[chan2]%4)];
            if(crosstalk[ge1][c2][bin] != -1 ){
              //  correction = crosstalk[ge1][c2][bin] + ((crosstalk[ge1][c2][bin+1] - crosstalk[ge1][c2][bin]) * (float)(((1940+dt)%160)/160));
              correction = crosstalk[ge1][c2][bin];
              //  fprintf(stdout,"CT enter, %d %d: %d %f %f %f %f: %f + %f = %f\n",chan,chan2,bin,crosstalk[ge1][c2][bin+1],crosstalk[ge1][c2][bin],(float)(((1940+dt)%160)/160),alt->ecal,ptr->ecal,(alt->ecal * correction),(ptr->ecal+(alt->ecal * correction)));
              ptr->ecal += alt->ecal * correction;
            }
          }
        }
      } // end of while

    } // end of if( ptr->subsys == SUBSYS_HPGE_A){
    return(0);
  }

// Presort - do Suppression and Addback here
//  - frag_idx is about to leave coinc window (which ends at end_idx)
//    check other frags in window for possible suppression and/or summing
//  also calculate multiplicities[store in frag_idx only]
int pre_sort_exit(int frag_idx, int end_idx)
{
  Grif_event *alt, *ptr = &grif_event[frag_idx];
  int i, dt;
  int chan,chan2;

  // Assign chan local variable and check it is a valid channel number
  if( (chan=ptr->chan)<0 || ptr->chan >= odb_daqsize ){
    fprintf(stderr,"presort error: ignored event in chan:%d\n",ptr->chan );
    return(-1);
  }
  i = frag_idx; ptr->multiplicity = 1;
  while( i != end_idx ){ // need at least two events in window
    if( ++i >=  PTR_BUFSIZE ){ i=0; } alt = &grif_event[i]; // WRAP
    if( (chan2=alt->chan)<0 || alt->chan >= odb_daqsize ){
      fprintf(stderr,"presort error: ignored event in chan:%d\n",alt->chan );
      continue;
    }

    // Determine absolute time difference between timestamps
    dt = ptr->ts - alt->ts; if( dt < 0 ){ dt = -1*dt; }

    // Restrict to 2 microseconds presort window for everything except HPGe crosstalk to maintain speed
    if(alt->subsys != SUBSYS_HPGE_A && dt>250){ continue; }

    // Determine multiplicity
    if( alt->subsys == ptr->subsys ){ ++ptr->multiplicity; }

    // SubSystem-specific pre-processing
    switch(ptr->subsys){
      case SUBSYS_HPGE_A:

      // HPGe pile-up corrections
      // THE PRE_SORT WINDOW SHOULD BE EXTENDED TO COVER THE FULL POSSIBLE TIME DIFFERENCE BETWEEN PILE-UP events
      // THIS IS EQUAL TO THE DIFF PERIOD OF HPGE TYPE
      // First assign the pileup class type, then correct the energies
      if(alt->subsys == SUBSYS_HPGE_A && chan2 == chan){
        perform_pileup_correction(ptr, alt, dt, chan, chan2, i, end_idx);
      }

      /*if(alt->subsys == SUBSYS_BGO){
        printf("crystals: [%i %i], dt: %i\n",crystal_table[ptr->chan],crystal_table[alt->chan],dt);
      }*/

      // BGO suppression of HPGe
      if( (dt >= bgo_window_min && dt <= bgo_window_max) && alt->subsys == SUBSYS_BGO && !ptr->suppress ){
        if( crystal_table[ptr->chan]/16 == crystal_table[alt->chan]/16 ){ ptr->suppress = 1; } //per-clover suppression
        //if( crystal_table[ptr->chan] == crystal_table[alt->chan] ){ ptr->suppress = 1; } //per-crystal suppression
      }
      /*// Germanium addback -
      //    earliest fragment has the sum energy, others are marked -1
      // Remember the other crystal channel number in alt_chan for use in Compton Polarimetry
      if( (dt >= addback_window_min && dt <= addback_window_max) && alt->subsys == SUBSYS_HPGE_A ){
        if( alt->esum >= 0 && crystal_table[alt->chan]/16 == crystal_table[ptr->chan]/16 ){
          ptr->esum += alt->esum; alt->esum = -1; ptr->alt_chan = alt->chan;
        }
      }*/
      break;
      case SUBSYS_HPGE_B:
      // HPGe B pile-up corrections
      if(alt->subsys == SUBSYS_HPGE_B && chan2 == chan){
        perform_pileup_correction(ptr, alt, dt, chan, chan2, i, end_idx);
      }
      break;
      default: break; // Unrecognized or unprocessed subsys type
    }// end of switch
  }// end of while

  return(0);
}

      // HPGe pile-up corrections
      // THE PRE_SORT WINDOW SHOULD BE EXTENDED TO COVER THE FULL POSSIBLE TIME DIFFERENCE BETWEEN PILE-UP events
      // THIS IS EQUAL TO THE DIFF PERIOD OF HPGE TYPE
      // First assign the pileup class type, then correct the energies
int perform_pileup_correction(Grif_event *ptr, Grif_event *alt, int dt, int chan, int chan2, int i, int end_idx)
{
  Grif_event *alt2;
  int j, dt13;
  float k1,k2, energy, correction, correction12, correction23;

  if(ptr->pileup==1 && ptr->nhit ==1){
    ptr->pu_class = PU_SINGLE_HIT; // no pileup, this is the most common type of HPGe event
    return(0);
  }else if(ptr->pileup==0){
    ptr->pu_class = PU_ERROR; // Pileup class, error
    return(0);
  }
if(dt>500){ return(0); } // Restrict pileup handling to 5 microseconds. 8 microseconds needed in some early datasets

  // Two hit pileup...
  if((ptr->pileup==1 && ptr->nhit==2) && (alt->pileup==2 && alt->nhit==1)){
    ptr->pu_class = alt->pu_class = PU_2HIT_ERROR; // Pileup class, error for 2Hits
    ptr->delta_t = alt->delta_t = dt;
    if(ptr->q1>0 && ptr->integ1>0 && ptr->q2>0 && ptr->integ2>0 && alt->q1>0 && alt->integ1>0){
      // 2 Hit, Type A. The (ptr) fragement is the first Hit of a two Hit pile-up event.
      ptr->pu_class = PU_2HIT_A1ST; alt->pu_class = PU_2HIT_A2ND; // Pileup class, 1st and 2nd of 2Hits
      ptr->delta_t = alt->delta_t = dt;  // time difference between hits
    }else{
      // 2 Hit, Type B, where 2nd Hit integration region starts after 1st Hit integration ends but before 2nd Hit CFD has completed
      ptr->pu_class = PU_2HIT_C1ST; alt->pu_class = PU_2HIT_C2ND; // Pileup class, 1st and 2nd of 2Hits
      ptr->delta_t = alt->delta_t = dt;  // time difference between hits
    }
  }else if((ptr->pileup==1 && ptr->nhit==2) && (alt->pileup==1 && alt->nhit==1)){
    // 2 Hit, Type C, where 2nd Hit integration region starts after 1st Hit integration and 2nd Hit CFD have ended
    ptr->pu_class = PU_2HIT_B1ST; alt->pu_class = PU_2HIT_B2ND; // Pileup class, 1st and 2nd  of 2Hits
    ptr->delta_t = alt->delta_t = dt; // Save the time difference between pileup hits into both hits
  }
 // Three hit pileup...
  else if((ptr->pileup==1 && ptr->nhit==3) && (alt->pileup==2 && alt->nhit==2)){ // 3Hit pileup
    ptr->pu_class = alt->pu_class = PU_3HIT_ERROR; // Pileup class, error for 3Hits
    if(ptr->q1>0 && ptr->integ1>0 && ptr->q2>0 && ptr->integ2>0 && alt->q1>1 && alt->integ1>0 && alt->q2>0 && alt->integ2>0){
      j=i+1;
      while( j != end_idx ){ // need to find the third events in window associated with this channel
        if( ++j >=  PTR_BUFSIZE ){ break; } alt2 = &grif_event[j]; // WRAP
        if(alt2->chan == chan){ // It must also be a HPGe if the channel number is the same
          if(alt2->pileup==3 && alt2->nhit==1){
            alt2->pu_class = PU_3HIT_3RD; // Pileup class
            if(alt2->q1>1 && alt2->integ1>0){
              // Determine absolute time difference between timestamps for Hit 1 and 3
              dt13 = ptr->ts - alt2->ts; if( dt13 < 0 ){ dt13 = -1*dt13; }
              // The Differencitation period of the HPGe Pulse Height evaluation is L = 5000ns.
              if(dt13>500){
                // Triple pileup case A ... in which the 3rd pulse occurs more than L samples after the first
                //                          there are 5 regions, 3 of which are not piled up (1 per pulse)
                correction23 = (alt->q1/alt->integ1)-((alt->q2/alt->integ2)-(alt2->q1/alt2->integ1));
                correction12 = (ptr->q1/ptr->integ1)-((ptr->q2/ptr->integ2)-(alt->q1/alt->integ1)-correction23);
                // Hit 1
                ptr->pu_class = PU_3HIT_1ST;
                energy = (spread(ptr->q1)/ptr->integ1) + correction12;
                ptr->ecal=ptr->esum = offsets[chan]+energy*(gains[chan]+energy*quads[chan]);
                // Hit 2
                alt->delta_t = dt; alt->pu_class = PU_3HIT_2ND;
                energy = (spread(alt->q1)/alt->integ1) - correction12 + correction23;
                alt->ecal=alt->esum = offsets[chan2]+energy*(gains[chan2]+energy*quads[chan2]);
                // Hit 3
                alt2->delta_t = dt13; alt2->pu_class = PU_3HIT_3RD;
                energy = (spread(alt2->q1)/alt2->integ1) - correction23;
                alt2->ecal=alt2->esum = offsets[alt2->chan]+energy*(gains[alt2->chan]+energy*quads[alt2->chan]);
              }else{
                // Triple pileup case B ... in which the 3rd pulse occurs less than L samples after the first
                //                          again 5 regions, only 2 of which are not piled up (first and last pulse)
                //                          There is no region to obtain the height of pulse 2
                //                          so the event contains K12, the sum of pulse 1+2, in place of pulseheight2
                correction23 = (alt->q1/alt->integ1)-((alt->q2/alt->integ2)-(alt2->q1/alt2->integ1));
                correction12 = (ptr->q1/ptr->integ1)-((ptr->q2/ptr->integ2)-(alt->q1/alt->integ1)-correction23);
                // Hit 1
                ptr->pu_class = PU_3HIT_1ST;
                energy = (spread(ptr->q1)/ptr->integ1) + correction12;
                ptr->ecal=ptr->esum = offsets[chan]+energy*(gains[chan]+energy*quads[chan]);
                // Hit 2
                alt->delta_t = dt; alt->pu_class = PU_3HIT_2ND;
                energy = (spread(alt->q1)/alt->integ1) - correction12 + correction23;
                alt->ecal=alt->esum = offsets[chan2]+energy*(gains[chan2]+energy*quads[chan2]);
                // Hit 3
                alt2->delta_t = dt13; alt2->pu_class = PU_3HIT_3RD;
                energy = (spread(alt2->q1)/alt2->integ1) - correction23;
                alt2->ecal=alt2->esum = offsets[alt2->chan]+energy*(gains[alt2->chan]+energy*quads[alt2->chan]);
              }
              break; // Break the while if we found the third Hit
            }
          }
        }
      } // end of while for triple coincidence
    }
  } // end of 3Hit pileup type assignments

  // Now apply hit-specific energy corrections
  if(pileupk1[chan][0] != 1){
    if(ptr->pu_class>=PU_2HIT_A1ST && ptr->pu_class<=PU_2HIT_C2ND){ // 2-Hit pileup events
      // Apply the k1 dependant correction to the energy of the first hit
      // It was already checked that chan for ptr and alt are the same for pileup events
      k1 = ptr->integ1;
      ptr->ecal=ptr->esum = ptr->ecal*( pileupk1[chan][0]+(k1*pileupk1[chan][1])+(k1*k1*pileupk1[chan][2])+(k1*k1*k1*pileupk1[chan][3])
      +(k1*k1*k1*k1*pileupk1[chan][4])+(k1*k1*k1*k1*k1*pileupk1[chan][5])+(k1*k1*k1*k1*k1*k1*pileupk1[chan][6]));

      // Apply the E1 and k2 dependant offset correction to the energy of the second hit
      // Apply the k2 dependant correction to the energy of the second hit
      k2 = alt->integ1;
      correction = ptr->ecal*( pileupE1[chan][0]+(k2*pileupE1[chan][1])+(k2*k2*pileupE1[chan][2])+(k2*k2*k2*pileupE1[chan][3])
      +(k2*k2*k2*k2*pileupE1[chan][4])+(k2*k2*k2*k2*k2*pileupE1[chan][5])+(k2*k2*k2*k2*k2*k2*pileupE1[chan][6]));
      alt->ecal=alt->esum = (alt->ecal*( pileupk2[chan][0]+(k2*pileupk2[chan][1])+(k2*k2*pileupk2[chan][2])+(k2*k2*k2*pileupk2[chan][3])
      +(k2*k2*k2*k2*pileupk2[chan][4])+(k2*k2*k2*k2*k2*pileupk2[chan][5])+(k2*k2*k2*k2*k2*k2*pileupk2[chan][6])))+correction;
    }
  }
  alt->alt_ecal=ptr->ecal; // Remember the ecal of the first Hit in this second Hit. Must be done regardless if a correction is made

  //printf("%f %f %f\n",pileupE1[chan][0],pileupE1[chan][1],pileupE1[chan][2]); //print some pileup correction coefficents

  return(0);
}

//writes data for a single sorted_evt in a SMOL tree
int lastWinIdx = -1;
uint64_t firstHitTs;
uint8_t fill_smol_entry(FILE *out, const int win_idx, const int frag_idx)
{
  //fprintf(stdout,"Called fill entry for win: %i, frag: %i, last win: %i\n",win_idx,frag_idx,lastWinIdx);
  if(win_idx == 0){
    lastWinIdx = -1; //wrap
  }
  if(frag_idx < win_idx){
    lastWinIdx = -1; //wrap
    return 0;
  }
  if(win_idx <= lastWinIdx){
    return 0; //don't double fill
  }
  Grif_event *ptr;
  int i;

  // initialize SMOL tree event
  sorted_evt *sortedEvt = calloc(1,sizeof(sorted_evt));
  sortedEvt->header.numHPGeHits = 0;

  for(i=win_idx; ; i++){
    
    ptr = &grif_event[i];

    if( i >= PTR_BUFSIZE ){ i=0; } //wrap 
    //if( i != win_idx && flag == SORT_ONE ){ break; }
    if( ptr->dtype == 15 ){ if( i==frag_idx ){ break; } continue; } // scalar
    if( ptr->chan == -1 ){
        printf("DefSort: UNKNOWN_CHAN type=%d\n", ptr->dtype);
        if( i==frag_idx ){ break; } continue;
    } //  ????
    
    /*if( ptr->chan<0 || ptr->chan >= odb_daqsize ){
      //fprintf(stderr,"SmolSort: UNKNOWN_CHAN=%i type=%d\n",ptr->chan,ptr->dtype);
      if( i==frag_idx ){ break; } continue;
    }
    if( i != win_idx && i==frag_idx ){ break; }
    if( ptr->dtype == 15 ){ if( i==frag_idx ){ break; } continue; } // scalar*/
    //fprintf(stdout,"  checking index %i (ts=%li)\n",i,ptr->ts);
    
    switch(ptr->subsys){
      case SUBSYS_HPGE_A: // Ge
        // Only use GRGa
        if(ptr->suppress != 1){
          //^passes Compton suppression
          if((ptr->pu_class >= 0)&&(ptr->pu_class < 16)){
            psd_vals[ptr->pu_class]++;
          }
          //if(ptr->pu_class == PU_SINGLE_HIT){
            //^no pileup (will want to include pileup correction later)
            int c1 = crystal_table[ptr->chan];
            if( c1 >= 0 && c1 < 64){
              if(sortedEvt->header.numHPGeHits >= MAX_EVT_HIT){
                fprintf(stderr,"WARNING: too many hits in win_idx %i, frag_idx %i\n",win_idx,frag_idx);
                break;
              }
              if((ptr->ecal > 5.0f)&&(ptr->ecal < 16384.0f)){ //try to filter out weird low and high energy stuff
                
                
                //Handle energy and position
                //printf("Energy: %f\n",ptr->ecal);
                sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].energy = ptr->ecal;
                sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].core = (uint8_t)(c1);
                if(sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].core >= 64){
                  fprintf(stderr,"WARNING: invalid GRIFFIN core: %u",sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].core);
                  break;
                }
                if(ptr->pu_class != PU_SINGLE_HIT){
                  sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].core |= (uint8_t)((uint8_t)(1) << 7); //set pileup flag
                }

                //Handle timing
                //has to be done AFTER position is set, since the CFD fail flag modifies the position value
                double grifT = getGrifTime(ptr);
                double CFDDiff = grifT - (double)(ptr->ts)*10.0; //difference between timestamp and CFD corrected time, in ns
                if(fabs(CFDDiff) > 300.0){
                  //printf("CFD fail, diff: %f\n",CFDDiff);
                  sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].core |= (uint8_t)((uint8_t)(1) << 6); //set CFD fail flag
                  grifT = (double)(ptr->ts)*10.0; //use timestamp time if the CFD fails
                }
                //printf("time: %f\n",grifT);

                if(grifT > 0.0){
                  if(sortedEvt->header.evtTimeNs == 0){
                    sortedEvt->header.evtTimeNs = grifT;
                    firstHitTs = ptr->ts;
                  }
                  sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].timeOffsetNs = (float)(grifT - sortedEvt->header.evtTimeNs);
                  
                  int tsDiff = (int)(ptr->ts - firstHitTs);
                  if(tsDiff < 0){
                    printf("Negative timestamp difference:\n",firstHitTs);
                    printf("  Hit 0: %lu\n",firstHitTs);
                    printf("  Hit %i: %lu\n",sortedEvt->header.numHPGeHits,ptr->ts);
                    sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].tsDiff = 255U;
                  }else if(tsDiff <= 255){
                    sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].tsDiff = (uint8_t)(tsDiff);
                  }else{
                    printf("Timestamp difference exceeding 255:\n",firstHitTs);
                    printf("  Hit 0: %lu\n",firstHitTs);
                    printf("  Hit %i: %lu\n",sortedEvt->header.numHPGeHits,ptr->ts);
                    sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].tsDiff = 255U;
                  }
                  //printf("tsDiff: %i\n",tsDiff);

                  //filter out duplicate data
                  uint8_t dupFound = 0;
                  for(int i = 0; i<sortedEvt->header.numHPGeHits;i++){
                    if((sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].core & 63U) == (sortedEvt->hpgeHit[i].core & 63U)){
                      if(fabs(sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].timeOffsetNs - sortedEvt->hpgeHit[i].timeOffsetNs) < 200.0){
                        dupFound = 1; //duplicate hit found
                      }
                    }
                    /*if(fabs(sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].timeOffsetNs - sortedEvt->hpgeHit[i].timeOffsetNs) < 1.0){
                      printf("  hit: %u, core: %u, time: %f ns, offset: %f ns\n",i,sortedEvt->hpgeHit[i].core,sortedEvt->header.evtTimeNs + sortedEvt->hpgeHit[i].timeOffsetNs,sortedEvt->hpgeHit[i].timeOffsetNs);
                      printf("  hit: %u, core: %u, time: %f ns, offset: %f ns\n",sortedEvt->header.numHPGeHits,sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].core,sortedEvt->header.evtTimeNs + sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].timeOffsetNs,sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].timeOffsetNs);
                    }*/
                  }
                  if(dupFound == 0){
                    //hit is not a duplicate, so store it
                    //printf("passed hit: %u, core: %u, energy: %0.3f, suppress: %i, ts: %li, cfd: %i, time: %0.2f ns\n",sortedEvt->header.numHPGeHits,sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].core,sortedEvt->hpgeHit[sortedEvt->header.numHPGeHits].energy,ptr->suppress,ptr->ts,ptr->cfd,grifT);
                    sortedEvt->header.numHPGeHits++;
                  }
                }
                
              }
            }else{
              fprintf(stderr,"WARNING: unknown GRIFFIN crystal %i\n",c1);
            }
          //}
        }
        break; // outer-switch-case-GE
      case SUBSYS_BGO:
        //at least one suppressor fired
        sortedEvt->header.metadata |= (uint8_t)(1U << 1);
        break;
      default:
        break; // Unrecognized or unprocessed subsys type
    }// end of switch(ptr)

    lastWinIdx = i;
    if( i==frag_idx ){ break; }
  }
  
  if((sortedEvt->header.numHPGeHits > 0)&&(sortedEvt->header.numHPGeHits <= MAX_EVT_HIT)){

    //finalize sorted event data
    sortedEvt->header.metadata |= (uint8_t)(1U << 7); //set data validation bit

    //write sorted event to SMOL tree
    fwrite(&sortedEvt->header,sizeof(evt_header),1,out);
    //write hits
    for(uint8_t j = 0; j<sortedEvt->header.numHPGeHits;j++){
      fwrite(&sortedEvt->hpgeHit[j].timeOffsetNs,sizeof(float),1,out);
      fwrite(&sortedEvt->hpgeHit[j].energy,sizeof(float),1,out);
      fwrite(&sortedEvt->hpgeHit[j].tsDiff,sizeof(uint8_t),1,out);
      fwrite(&sortedEvt->hpgeHit[j].core,sizeof(uint8_t),1,out);
      //fprintf(stdout,"Hit %u - core: %u, energy: %0.2f, time offset: %0.2f, win: %i, frag: %i\n",j,sortedEvt->hpgeHit[j].core,(double)sortedEvt->hpgeHit[j].energy,(double)sortedEvt->hpgeHit[j].timeOffsetNs,win_idx,frag_idx);
    }
    free(sortedEvt);
    return 1; //sorted an event
  }

  free(sortedEvt);
  return 0;
  
}

//#######################################################################
//###########   READ XML ODB DUMP FROM START OF DATA FILE   #############
//#######################################################################

// Note - the odb names do not distinguish between subtypes of detectors
// e.g Ge A and B channels
// the subsystem names will be extended to include this information
// (and the odb-specific names are only used below)

#define MAX_ODB_SUBSYS 24
#define ODBHANDLE_GRG   0
#define ODBHANDLE_GRS   1
#define ODBHANDLE_SEP   2
#define ODBHANDLE_PAC   3
#define ODBHANDLE_LBS   4
#define ODBHANDLE_LBT   5
#define ODBHANDLE_LBL   6
#define ODBHANDLE_DSC   7
#define ODBHANDLE_ART   8
#define ODBHANDLE_ZDS   9
#define ODBHANDLE_RCS  10
#define ODBHANDLE_XXX  11
#define ODBHANDLE_DSW  12
#define ODBHANDLE_DSG  13
#define ODBHANDLE_DAL  14
#define ODBHANDLE_DAT  15
#define ODBHANDLE_QED  16
#define ODBHANDLE_UNK  23
static char odb_handle[MAX_ODB_SUBSYS][8] = {
  "GRG", "GRS", "SEP", "PAC",  //  0- 3
  "LBS", "LBT", "LBL", "DSC",  //  4- 7
  "ART", "ZDS", "RCS", "XXX",  //  8-11
  "DSW", "DSG", "DAL", "DAT",  //  12-15
  "QED",    "",    "",    "",
  "",    "",    "",    "UNK"
};

static char   path[256];
static char dirname[64],value[32];
extern char midas_runtitle[SYS_PATH_LENGTH];

static void *arrayptr;
int read_odb_items(int len, int *bank_data)
{
   char *path_ptr, *ptr, *str, *odb_data = (char *)bank_data;
   int i=0, c = '<', d = '>', dtype=0, active=0, index=0;
   ptr = odb_data;  path_ptr = path;
   while(1){
      if( (str = strchr(ptr,c)) == NULL ){ break; }
      ptr = str;
      if( (str = strchr(ptr,d)) == NULL ){ break; }

      if( strncmp(ptr,"<!--",4) == 0 || strncmp(ptr,"<odb", 4) == 0 ||
                                        strncmp(ptr,"</odb",5) == 0 ){ // comment - skip
      } else if( strncmp(ptr,"<dir ",5) == 0 ){
         if( strncmp(ptr,"<dir name=\"",11) == 0 ){
            i=11; while( *(ptr+i) != '"' && *(ptr+i) != d ){ ++i; }
         }
         memcpy(dirname, ptr+11, i-11); dirname[i-11] = '\0';
         if( *(ptr+1+i) == '/' ){ ptr=str+1; continue; }
         //if( sscanf(ptr,"<dir name=\"%s\">", dirname) < 1 ){
         //   fprintf(stderr,"can't read dirname\n"); ptr=str+1; continue;
         //}
         //if( strncmp(dirname+strlen(dirname)-3,"\"/>",3) == 0 ){
         //   ptr=str+1; continue;
         //}
         //if( dirname[strlen(dirname)-1]=='>'  ){
         //   dirname[strlen(dirname)-1]='\0';
         //}
         //if( dirname[strlen(dirname)-1]=='\"' ){
         //  dirname[strlen(dirname)-1]='\0';
         //}
         *path_ptr = '/'; strcpy(path_ptr+1, dirname);
         path_ptr += strlen(dirname)+1;
         *path_ptr = '\0';
      } else if( strncmp(ptr,"</dir>",6) == 0 ){
         while(1){
            if( --path_ptr < path ){ path_ptr = path;  *path_ptr = '\0';  break; }
            if( *path_ptr == '/' ){ *path_ptr = '\0';  break; }
         }
         index=0; // for debugger to stop here
      } else if( strncasecmp(ptr,"<key name=\"Run Title\" type=\"STRING\"", 35) == 0 ){
         ptr = str+1;
         if( (str = strchr(ptr,c)) == NULL ){ break; }
         i = (str-ptr) > SYS_PATH_LENGTH-1 ? SYS_PATH_LENGTH-1 : (str-ptr);
         memcpy( midas_runtitle, ptr, i ); midas_runtitle[i] = 0;
         ptr += i+1;
         if( (str = strchr(ptr,d)) == NULL ){ break; }
      } else if( strncmp(ptr,"</keyarray>",10) == 0 ){ active = 0; arrayptr = (void *)('\0');
      } else if( strncmp(ptr,"<keyarray ",10) == 0 ){
         if( strcmp(path,"/DAQ/params/MSC") != 0 &&
             strcmp(path,"/DAQ/MSC")        != 0 &&
             strcmp(path,"/DAQ/PSC")        != 0 ){  ptr=str+1; continue; }
         if( sscanf(ptr,"<keyarray name=\"%s", value) < 1 ){
            fprintf(stderr,"can't read keyarray entry\n"); ptr=str+1; continue;
         }
         if( value[strlen(value)-1]=='\"' ){ value[strlen(value)-1]='\0'; }
         if( strcmp(value,"PSC") == 0 || strcmp(value,"MSC") == 0 ){
            active = 1; arrayptr = (void *)addr_table; dtype=1;
         }
         if( strcmp(value,"chan") == 0 ){
            active = 1; arrayptr = (void *)chan_name; dtype=3;
         }
         if( strcmp(value,"datatype") == 0 ){
          //active = 1; arrayptr = (void *)dtype_table; dtype=1;
            active = 1; arrayptr = (void *)dtype_table; dtype=0;
         }
         if( strcmp(value,"gain") == 0 ){
            active = 1; arrayptr = (void *)gain_table; dtype=2;
         }
         if( strcmp(value,"offset") == 0 ){
            active = 1; arrayptr = (void *)offs_table; dtype=2;
         }
         if( strcmp(value,"quadratic") == 0 ){
            active = 1; arrayptr = (void *)quad_table; dtype=2;
         }
      } else if( strncmp(ptr,"<value index=",13) == 0 ){
         if( !active ){ ptr=str+1; continue; }
         // remove the >< surrounding the value, and move str to the end of the line
         *str = ' '; if( (str = strchr(str,c)) == NULL ){ break; }
         *str = ' '; if( (str = strchr(str,d)) == NULL ){ break; }
         if( sscanf(ptr,"<value index=\"%d\" %s /value>", &index, value) < 2 ){
            fprintf(stderr,"can't read value entry\n");
         }
         if( index < 0 || index >= MAX_DAQSIZE ){
            fprintf(stderr,"index %d out of range\n", index);
         }
         // index starts at zero, odb_daqsize is count
         if( index >= odb_daqsize ){ odb_daqsize = index+1; }
         if(        dtype == 0 ){  // int
            if( sscanf(value,"%d", (((int *)arrayptr)+index)) < 1 ){
               fprintf(stderr,"can't read value %s\n", value);
            }
         } else if( dtype == 1 ){  // short int
            if( sscanf(value,"%hd", (((short *)arrayptr)+index)) < 1 ){
               fprintf(stderr,"can't read value %s\n", value);
            }
         } else if( dtype == 2 ){  // float
            if( sscanf(value,"%f", (((float *)arrayptr)+index)) < 1 ){
               fprintf(stderr,"can't read value %s\n", value);
            }
         } else {                 // string
            strncpy(arrayptr+index*CHAN_NAMELEN, value, CHAN_NAMELEN);
            *((char *)arrayptr+(index+1)*CHAN_NAMELEN - 1) = '\0';
         }
      }
      ptr=str+1;
   }
   fprintf(stdout,"odb record: %d bytes\n", len);

   // arrays typically around 500 entries [one per "chan"] each entry with ...
   //   daq-address, name, type, gains etc.
   //
   gen_derived_odb_tables();

   return(0);
}

// original odb arrays were read into {addr_table,chan_name,dtype_table(+gains)}
// extract extra details stored in channel names (and record for later)
// (these details include crystal/element numbers and polarities)
// [use above for subsystem (no longer use datatype to determine subsystems)]
extern int read_caen_odb_addresses(int odb_daqsize, unsigned short *addr_table);
int gen_derived_odb_tables()
{
   int i, j, tmp, subsys, pos, element, output_type;
   char sys_name[64], crystal, polarity, type;

   read_caen_odb_addresses(odb_daqsize, (unsigned short *)addrs);

   // generate reverse mapping of address to channel number
   //  (most of this array is undefined and stays at -1)
   memset(address_chan, 0xFF, sizeof(address_chan)); // set to -1
   for(i=0; i<MAX_ADDRESS && i<odb_daqsize; i++){
      address_chan[ (unsigned short)chan_address[i] ] = i;
   }

   memset(crystal_table,  0xff, MAX_DAQSIZE*sizeof(int)); // initialise all to -1
   memset(element_table,  0xff, MAX_DAQSIZE*sizeof(int));
   memset(polarity_table, 0xff, MAX_DAQSIZE*sizeof(int));
   memset(subsys_table,   0xff, MAX_DAQSIZE*sizeof(int));
   for(i=0; i<MAX_DAQSIZE && i<odb_daqsize; i++){
      if( (tmp=sscanf(chan_name[i], "%3c%d%c%c%d%c", sys_name, &pos, &crystal, &polarity, &element, &type)) != 6 ){
         fprintf(stderr,"can't decode name[%s] decoded %d of 6 items\n", chan_name[i], tmp );
         continue;
      }
      for(j=0; j<MAX_ODB_SUBSYS; j++){
         if( strncmp(sys_name, odb_handle[j], 3) == 0 ){ subsys = j; break; }
      }
      if( j == MAX_ODB_SUBSYS ){ subsys = j-1; // use final entry: "unknown"
         fprintf(stderr,"Unknown subsystem[%s] in %s\n", sys_name, chan_name[i]);
      }

      // Mention bad detector types (no longer relied on for subsystem id)
      if( dtype_table[i] < 0 || dtype_table[i] >= 16 ){
         fprintf(stderr,"bad datatype[%d] at table position %d\n", dtype_table[i], i);
      }

      // Some detector elements have more than one output (HPGe A and B)
      // 1 is A, 0 is B, -1 is X or unknown
      output_type = type=='A' ? 1 : (type=='B' ? 0 : -1);

      // Polarity: 1 is N, 0 is P or T or S, -1 is anything else
      if(        polarity == 'N' ){ polarity_table[i] = 1;
      } else if( polarity == 'P' ){ polarity_table[i] = 0;
      } else if( polarity == 'T' ){ polarity_table[i] = 0; // TAC signal
      } else if( polarity == 'S' ){ polarity_table[i] = 1; // ARIES Standard Ouput signal
      } else if( polarity == 'F' ){ polarity_table[i] = 0; // ARIES Fast Output signal
      } else if( polarity == 'X' ){ polarity_table[i] = 0; // XXX type
      } else { fprintf(stderr,"unknown polarity[=%c] in %s\n", polarity, chan_name[i]); }

      // Record crystal and element numbers [** Naming schemes are subsystem-dependant **]
      switch(subsys){
      case ODBHANDLE_LBL: case ODBHANDLE_LBS: // LaBr,Paces, Ares and Zds
      case ODBHANDLE_LBT: case ODBHANDLE_ART:
      case ODBHANDLE_DAL: case ODBHANDLE_DAT:
      case ODBHANDLE_PAC: case ODBHANDLE_ZDS: case ODBHANDLE_DSW:
         crystal_table[i] = pos;
         if(        crystal == 'A' ){ element_table[i] = 1;
         } else if( crystal == 'B' ){ element_table[i] = 2;
         } else if( crystal == 'C' ){ element_table[i] = 3;
         } else if( crystal == 'X' ){ element_table[i] = -1; // just one crystal for LaBr3, ZDS, ART
         } else {
            fprintf(stderr,"unknown crystal for ancillary[=%c] in %s\n", crystal, chan_name[i]);
         } break;
      case ODBHANDLE_RCS:
         crystal_table[i] = pos;
         element_table[i] = element;
         break;
      case ODBHANDLE_GRG: case ODBHANDLE_GRS:
         element_table[i] = element;
         pos -= 1; pos *=4;
         if(        crystal == 'B' ){ crystal_table[i] = pos;
         } else if( crystal == 'G' ){ crystal_table[i] = pos+1;
         } else if( crystal == 'R' ){ crystal_table[i] = pos+2;
         } else if( crystal == 'W' ){ crystal_table[i] = pos+3;
         } else if( crystal == 'X' ){ crystal_table[i] = -1; // crystal undefined
         } else {
            fprintf(stderr,"unknown crystal[=%c] in %s\n", crystal, chan_name[i]);
         } break;
      default: break;
      }

      //printf("Channel %i name: %s - element: %i - crystal: %i\n",i,chan_name[i],element_table[i],crystal_table[i]); //debug

      // set full subsystem id (including polarity/output-type etc)
      switch(subsys){
      case ODBHANDLE_GRS: subsys_table[i] = SUBSYS_BGO;       break;
      case ODBHANDLE_SEP: subsys_table[i] = SUBSYS_SCEPTAR;   break;
      case ODBHANDLE_PAC: subsys_table[i] = SUBSYS_PACES;     break;
      case ODBHANDLE_LBS: subsys_table[i] = SUBSYS_LABR_BGO;  break;
      case ODBHANDLE_LBL: subsys_table[i] = SUBSYS_LABR_L;    break;
      case ODBHANDLE_DAL: subsys_table[i] = SUBSYS_LABR_L;    break;
      case ODBHANDLE_DSC: subsys_table[i] = SUBSYS_DESCANT;   break;
      case ODBHANDLE_RCS: subsys_table[i] = SUBSYS_RCMP;      break;
      case ODBHANDLE_DSW: subsys_table[i] = SUBSYS_DESWALL;  break;
      case ODBHANDLE_DSG: subsys_table[i] = SUBSYS_DSG;  break;
      case ODBHANDLE_GRG: subsys_table[i] = (output_type == 1) ? SUBSYS_HPGE_A :SUBSYS_HPGE_B; break;
      case ODBHANDLE_ZDS: subsys_table[i] = (output_type == 1) ? SUBSYS_ZDS_A  :SUBSYS_ZDS_B;  break;
      case ODBHANDLE_ART: subsys_table[i] = (polarity_table[i] == 1) ? SUBSYS_ARIES_A:SUBSYS_ARIES_B;break;
      case ODBHANDLE_XXX: subsys_table[i] = SUBSYS_IGNORE;    break;
      case ODBHANDLE_UNK: subsys_table[i] = SUBSYS_UNKNOWN;   break;
      case ODBHANDLE_DAT: if(crystal_table[i]<8){ subsys_table[i] = SUBSYS_TAC_LABR;
                          }else{ subsys_table[i] = SUBSYS_TAC_ZDS; }
                          break;
      case ODBHANDLE_LBT: if(crystal_table[i]<8){ subsys_table[i] = SUBSYS_TAC_LABR;
                          }else if(crystal_table[i]>8){ subsys_table[i] = SUBSYS_TAC_ART;
                          }else{ subsys_table[i] = SUBSYS_TAC_ZDS; }
                          break;
      }
   }
   memset(subsys_initialized, 0, sizeof(int)*MAX_SUBSYS );

   return(0);
}