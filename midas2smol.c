/* read midas data files
 set sources=( midas2smol.c midas-format.c grif-format.c config.c reorder.c user_sort.c default_sort.c )
 gcc -g     -o midas2smol $sources -rdynamic -ldl -lm -lpthread
 gcc -g -O3 -o midas2smol $sources -rdynamic -ldl -lm -lpthread
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include "config.h"
#include "grif-format.h"
#include "midas-format.h"

int sort_next_file(Config *cfg, Sort_status *sort);

int presort_events_cutoff = 256; // presort window typically much larger, and less processing is done
int coinc_events_cutoff = 256; //64 isn't enough for some high rate data, end up with high WinFull %s...
char midas_runtitle[SYS_PATH_LENGTH];
Sort_metrics diagnostics;
static Sort_status sort_status;
FILE *output_tree;
uint64_t psd_vals[MAX_PSD_VALS];
int main(int argc, char *argv[])
{
   if(argc < 3){
      fprintf(stdout,"  midas2smol midas_file output_SMOL_file config_json\n\n");
      fprintf(stdout,"    Sorts MIDAS files into the SMOL event list format.\n\n");
      fprintf(stdout,"    MIDAS files are expected to be named using the standard run\n");
      fprintf(stdout,"    and subrun numbering scheme (eg. run29623_000.mid).\n");
      fprintf(stdout,"    Provide the path of the first subrun, and all available\n");
      fprintf(stdout,"    subruns for that run number will then be sorted.\n\n");
      fprintf(stdout,"    If no config file name is given, the default name\n");
      fprintf(stdout,"    (last.json) will be assumed.\n\n");
      exit(0);
   }
   
   srand(28719747); //seed random number generator (use fixed seed so that results are consistent between sorts)
   Sort_status *sort = &sort_status;
   Config *cfg;
   if(argc >= 4){
      init_config(argv[3]);
   }else{
      init_config(DEFAULT_CONFIG);
   }
   copy_config(configs[0], configs[1]); // copy config0 to cfg1 for sorting
   cfg = configs[1];
   sort->subrun = 0;
   sort->single_thread = 0;

   strncpy(cfg->out_file,argv[2],SYS_PATH_LENGTH-1); //setup output filename
   add_sortfile(argv[1]);
   if(sort->num_subruns > 1){
      fprintf(stdout,"MIDAS file: %s (%i subruns)\n",argv[1],sort->num_subruns);
   }else{
      fprintf(stdout,"MIDAS file: %s\n",argv[1]);
   }
   fprintf(stdout,"Output file: %s\n",cfg->out_file);
   
   if( open_next_sortfiles(sort) == 0 ){
      sort_next_file(cfg, sort);
      fprintf(stdout,"DONE\n");
      close_sortfiles(sort);
   }
   
   exit(0);
}

static int presort_window_start, sort_window_start;
static int done_events;
extern void grif_main(Sort_status *arg);
extern void reorder_main(Sort_status *arg);
extern void reorder_out(Sort_status *arg);
extern void sort_main(Sort_status *arg, FILE *out);
static pthread_t midas_thread, grif_thread, ordthrd, ordthr2;
static int reorder_save, singlethread_save, sortthread_save;
uint64_t numSortedEvts;
int sort_next_file(Config *cfg, Sort_status *sort)
{
   FILE *smolfp;
   time_t end, start=time(NULL);
   done_events = 0;
   presort_window_start = sort_window_start = 0;
   memset(&diagnostics, 0, sizeof(Sort_metrics) );
   diagnostics.run_sort_start = start;
   sort->shutdown_midas = sort->end_of_data = 0;
   sort->reorder_in_done = sort->reorder_out_done = 0;
   sort->grif_sort_done = sort->odb_done = sort->odb_ready = 0;
   reorder_save = sort->reorder;
   singlethread_save = sort->single_thread;
   sortthread_save = sort->sort_thread;
   numSortedEvts = 0UL;
   if( singlethread_save == 1 ){
      midas_main(sort);
   } else {
      printf("creating midas thread\n");
      pthread_create(&midas_thread, NULL, (void* (*)(void*))midas_main, sort);
      printf("creating reorder threads\n");
      pthread_create(&ordthrd,NULL,(void* (*)(void*))reorder_main, sort);
      pthread_create(&ordthr2,NULL,(void* (*)(void*))reorder_out, sort);

      while( !sort->odb_ready ){     // wait for midas thread to read odb event
         usleep(1);
      }
      
      init_default_sort(configs[1], sort);     // depend on odb in datafile
      pthread_create(&grif_thread, NULL,(void* (*)(void*))grif_main, sort);

      if( (smolfp=fopen(cfg->out_file,"wb")) != NULL ){
         printf("Opened output file: %s\n",cfg->out_file);
         memset(psd_vals,0,sizeof(psd_vals));
         uint64_t tmp = 0U; //placeholder
         fwrite(&tmp,sizeof(uint64_t),1,smolfp);
         fwrite(&psd_vals,sizeof(psd_vals),1,smolfp);
         sort_main(sort,smolfp); // this exits when sort is done
         //number of sorted events in SMOL format can only be 48 bits
         if(numSortedEvts > 0xFFFFFFFFFFFF){
            printf("WARNING: number of output events (%10ld) truncated to %10ld.\n",numSortedEvts,(uint64_t)(numSortedEvts & 0xFFFFFFFFFFFF));
         }
         numSortedEvts &= 0xFFFFFFFFFFFF;
         uint64_t smolFormatVersion = 1;
         numSortedEvts |= (smolFormatVersion << 48);
         //rebuild header
         fseek(smolfp,0,SEEK_SET);
         fwrite(&numSortedEvts,sizeof(uint64_t),1,smolfp);
         fwrite(&psd_vals,sizeof(psd_vals),1,smolfp);
         fclose(smolfp);
         printf("Wrote %10lu separated events to output file: %s\n",numSortedEvts & 0xFFFFFFFFFFFF,cfg->out_file);
      } else {
         printf("Can't open SMOL tree: %s to write\n",cfg->out_file);
         return(0);
      }

      sort->shutdown_midas = 1;
      pthread_join(midas_thread, NULL);
      pthread_join(ordthrd, NULL);
      pthread_join(ordthr2, NULL);
      pthread_join(grif_thread, NULL);
   }
   end=time(NULL);
   cfg->midas_start_time = diagnostics.midas_run_start;
   cfg->midas_runtime    = diagnostics.midas_last_timestamp+1;
   cfg->midas_runtime   -= cfg->midas_start_time;
   memcpy(cfg->midas_title, midas_runtitle, SYS_PATH_LENGTH);
   printf("Sorting took %ld seconds (%.3f hours)\n", end-start, (end-start)/3600.0);
   return 0;
}

extern volatile long grifevent_wrpos;
volatile long grifevent_rdpos;
extern Grif_event grif_event[PTR_BUFSIZE];

extern volatile unsigned long bankbuf_wrpos;
extern volatile unsigned long bankbuf_rdpos;
extern volatile long tsevents_in;
extern long tsevents_out;
extern volatile unsigned long eventbuf_rdpos;
extern volatile unsigned long eventbuf_wrpos;
extern volatile long grif_evcount;
void show_sort_state()
{
   int val = sort_status.midas_bytes/1000;
   int v2 = bankbuf_wrpos/1000000, v3 = bankbuf_rdpos/1000000;
   int v4 = tsevents_in/1000, v5 = tsevents_out/1000;
   int v6 = eventbuf_wrpos/1000000, v7 = eventbuf_rdpos/1000000;
   int v8 = grifevent_wrpos/1000, v9 = grifevent_rdpos/1000;
   int v10 = grifevent_wrpos - grifevent_rdpos;
   printf("MIDAS:read %d Mbytes [~%d Kevents]\n", val/1000, val/50);
   printf("      BUF in:%dMbytes out:%dMbytes  [Cap:%5.1f%%]\n",
          4*v2, 4*v3, (100000000.0*(v2-v3))/BANK_BUFSIZE );
   printf("REORDER: In:%dKevents Out%dKevents\n", v4, v5 );
   printf("     BUF In:%dMbytes  Out:%dMbytes  [Cap:%5.1f%%]\n",
          4*v6, 4*v7, (100000000.0*(v6-v7))/EVENTBUFSIZE );
   printf("GRIF: Unpacked:%dKevents  ", (int)(grif_evcount/1000) );
   printf("      Sorted  :%dKevents\n", done_events/1000 );
   printf("     BUF In:%dKevents Out%dKevents [=%d][Cap:%5.1f%%]\n\n",
          v8, v9, v10, (100.0*v10)/PTR_BUFSIZE);
}
Sort_status *get_sort_status(){ return( &sort_status ); }

void sort_main(Sort_status *arg, FILE *out)
{
   int nxtpos, rd_avail;
   static long grifevent_nxtpos;
   unsigned int usecs=100;

   printf("starting sort_main ...\n");
   grifevent_rdpos = grifevent_nxtpos = nxtpos = 0;
   while(1){
      // if( arg->shutdown_midas != 0 ){  break; }
      rd_avail = grifevent_wrpos - grifevent_nxtpos;
      if( arg->grif_sort_done && rd_avail < 1 ){ break; }
      if( rd_avail < 1 ){ usleep(usecs); continue; }
      numSortedEvts += process_event(&grif_event[nxtpos], nxtpos, out);
      //printf("sorted events: %lu\n",numSortedEvts);
      nxtpos = ++grifevent_nxtpos % PTR_BUFSIZE;
   }
   printf("sort_main finished\n");
}

static long calls, sorted, skipped, prefull, sortfull, completed_events;
// called when each new event read into ptr -> list[slot]
uint64_t process_event(Grif_event *ptr, int slot, FILE *out)
{
   time_t cur_time = time(NULL);
   static long prv_call;
   static time_t prv_time;
   int dt = cur_time - prv_time, de = calls - prv_call;
   if( prv_time == 0 ){ prv_time = cur_time; }
   calls++;

   if( cur_time-prv_time >= 10 ){
      prv_call = calls;
      printf("----------------------------------------------------------------\n");
      midas_status(cur_time); reorder_status(cur_time);  grif_status(cur_time);
      printf("ProcEvt: %10ld[Good:%3d%% Skip:%3d%% WinFull:%3d%%] %6.3f Mevt/s\n",
             calls, (int)(100.0*(calls-skipped-prefull)/calls),
             (int)(100.0*skipped/calls), (int)(100.0*prefull/calls), (de/(1000000.0*dt))
      );
      printf("SortedEvt: %10ld, CompletedEvt: %10ld\n",numSortedEvts,sorted);
      prv_time = cur_time;
   }

   return insert_presort_win(ptr, slot, out);
}

char *debug_show_ts(long ts)
{
   static char tmp[32];
   int deci_ms = ts/10000;
   int remain = ts - (deci_ms*10000);

   sprintf(tmp,"%5.1fms+%04d", 00000.1*deci_ms, remain);
   return(tmp);
}
char *debug_show_chan(Grif_event *ptr)
{
   static char tmp[32];
   if( ptr->chan != -1 ){
      sprintf(tmp," %3d", ptr->chan);
   } else {
      sprintf(tmp,"%04x", ptr->address);
   }
   return(tmp);
}
extern char chan_name[MAX_DAQSIZE][CHAN_NAMELEN];

// add event to presort window (event has just been read in)
//    recalculate coincwin (sorting any events leaving window)
// => final win of run won't be sorted, as these events will not leave window
uint64_t insert_presort_win(Grif_event *ptr, int slot, FILE *out)
{
   //printf("start insert_presort_win: %i %i\n",presort_window_start, slot);
   uint64_t numSort = 0;
   int win_count, win_end;
   Grif_event *alt;
   long dt;

   ///////////////// Presort window (used for suppression/addback)
   while( presort_window_start != slot ){
      //printf("while loop: %i %i\n",presort_window_start, slot);
      alt = &grif_event[presort_window_start];
      win_count = (slot - presort_window_start+2*PTR_BUFSIZE) % PTR_BUFSIZE;
      dt = ptr->ts - alt->ts; if( dt < 0 ){ dt *= -1; }

      // should exit while-loop when no more events outside window
      //    *BUT* add error recovery - if window too full, dump events
      if( dt < presort_window_width ){
         if( win_count < presort_events_cutoff ){ break; } // LIMIT to ?? events
         else { ++prefull; }
      }

      // event[win_start] is leaving window
      //    ( either because dt > coincwidth OR due to error recovery)
      // NOTE event[slot] is out of window - use slot-1 as window-end
      if( (win_end = slot-1) < 0 ){ win_end = PTR_BUFSIZE-1; } // WRAP
      //printf("start pre_sort_exit: %i %i\n",presort_window_start, win_end);
      pre_sort_exit(presort_window_start, win_end);
      //printf("end pre_sort_exit\n");
      numSort += insert_sort_win(alt, presort_window_start, out); // add event to next window
      if( ++presort_window_start >= PTR_BUFSIZE ){ presort_window_start=0; } // WRAP
      //printf("end while loop iter: %i %i %i\n",presort_window_start, slot, PTR_BUFSIZE);
   }
   // all events outside window have now been removed ...
   pre_sort_enter(presort_window_start, slot);
   //printf("numSort: %lu\n", numSort);
   return numSort;
}

// add event to main sort window (event has just left presort window)
uint64_t insert_sort_win(Grif_event *ptr, int slot, FILE *out)
{
   
   uint64_t numSort = 0;
   int win_count, win_end;
   Grif_event *alt;
   long dt;

   while( sort_window_start != slot ){ alt = &grif_event[sort_window_start];
       win_count = (slot - sort_window_start+2*PTR_BUFSIZE) % PTR_BUFSIZE;
      dt = ptr->ts - alt->ts; if( dt < 0 ){ dt *= -1; }

      // should exit while-loop when no more events outside window
      //    *BUT* add error recovery - if window too full, dump events
      if( dt < sort_window_width ){
       //if( win_count >= 0.45*PTR_BUFSIZE ){ ++sortfull; } else {
         if( win_count > coinc_events_cutoff ){ ++sortfull; } else {
            break;
         }
      }

      // event[win_start] is leaving window
      //    ( either because dt > coincwidth OR due to error recovery)
      // NOTE event[slot] is out of window - use slot-1 as window-end
      if( (win_end = slot-1) < 0 ){ win_end = PTR_BUFSIZE-1; } // WRAP
      if( alt->chan != -1 ){ ++sorted;
         numSort += (uint64_t)fill_smol_entry(out, sort_window_start, win_end);
      } else { ++skipped; }
      if( ++sort_window_start >= PTR_BUFSIZE ){ sort_window_start=0; } // WRAP
      ++grifevent_rdpos;  ++completed_events;
   }
   return numSort;
}
