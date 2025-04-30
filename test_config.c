// this is mostly pointless - just save everything into a "InitialConfig.cfg"
// (which is loaded on startup) and get rid of test_config.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grif-format.h"
#include "config.h"

//char *initvarlist[2*NUM_INITVARS]={
//   "GeE",       "Singles Ge Energy",               //  0
//   "GeSupE",    "Suppressed Singles Ge Energy",    //  1
//   "AbGeE",     "Addback Ge Energy",               //  2
//   "AbGeSupE",  "Suppressed Addback Ge Energy",    //  3
//   "GeClov",    "Ge Clover Number",                //  4
//   "GeCrystal", "Ge Crystal Number",               //  5
//   "GeT",       "Ge Cfd Time",                     //  6
//   "GeTs",      "Ge Timestamp",                    //  7
//   "SepE",      "Sceptar Energy",                  //  8
//   "SepT",      "Sceptar Time",                    //  9
//   "GeSepdT",   "Ge Sceptar time difference",      // 10
//   "GeSepdTs",  "Ge Sceptar timestamp difference", // 11
//   "GeGedT",    "Ge Ge time difference",           // 12
//};


#define NUM_INITGLOBALS 3
char *initgloballist[3*NUM_INITGLOBALS]={
   "Gamma-BGO Time difference",     "-50",   "50",
   "Gamma-Gamma Time difference",  "-100",  "100",
   "Beta-Gamma Time difference",    "-80",   "80",
};

#define NUM_INITGATES 3
char *initgatelist[4*NUM_INITGATES]={
   "BetaVeto",         "SEPE",  "=",   "0",
   "GeBetaDt",   "TD_GRG_SEP",  "<", "100",
   "BetaCoinc",        "SEPE",  ">",   "0",
};

#define NUM_INITHISTOS 3
char *inithistolist[5*NUM_INITHISTOS]={
   "GeEx", "GeBetaVeto",        "User/Test", "256", "HPGeEU",
   "GeE",  "Singles Ge Energy", "User/Test", "256", "HPGeEU",
   "SepT", "Sceptar Time",      "User/Test", "256", "SEPT",
};

#define NUM_INITGROUPS 3
char *initgrouplist[4*NUM_INITGROUPS]={
   "BetaVeto",      "BetaVeto",    "",          "",
   "GeBetaDt",      "BetaCoinc",   "GeBetaDt",  "",
   "BetaCoinc",     "BetaCoinc",   "",          "",
};

int init_default_config(Config *cfg) // setup initial live config (configs[0])
{
   static int done;
   int i, j, value;
   char *tmp;
   
   if( !done ){ done = 1; } else { return(0); } // only do once
   cfg->lock = 1;
   for(i=0; i< NUM_INITGLOBALS; i++){
      tmp = initgloballist[3*i+1];
      if( sscanf(tmp, "%d", &value) < 1 ){
         fprintf(stderr,"init_globals: can't read value in %s\n", tmp);
         cfg->lock = 0; return(-1);
      }
      tmp = initgloballist[3*i+2];
      if( sscanf(tmp, "%d", &j) < 1 ){
         fprintf(stderr,"init_globals: can't read value in %s\n", tmp);
         cfg->lock = 0; return(-1);
      }
      if(add_global(cfg,initgloballist[3*i],value,j)){cfg->lock=0;return(-1);}
   }
   cfg->lock = 0; return(0);
}
