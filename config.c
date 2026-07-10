#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <ctype.h>
#include "config.h"
#include "midas2smol.h"

///////////////////////////////////////////////////////////////////////////
/////////////////            Config  I/O            ///////////////////////
///////////////////////////////////////////////////////////////////////////
Config *configs[MAX_CONFIGS];

static char config_data[1024*1024];
int load_config(Config *cfg, const char *filename, char *buffer)
{
   int i,j, len, value, val2, val3, val4, val5, val6, address, type, instring;
  char *ptr, *name, *valstr, *title, *path, *var, *var2, op[8], tmp[80];
  float gain, offset, quad;
  float puk1[7], puk2[7], puE1[7];
  float ct0[16], ct1[16], ct2[16];
  // Initialize values to defaults
  // Values of -1 are ignored by edit_calibration - use this for all channels that are not HPGe to avoid bloating the size of the config
  float puk_reset[7]={1,0,0,0,0,0,0}, puE1_reset[7]={0,0,0,0,0,0,0}, pu_ignore[7]={-1,-1,-1,-1,-1,-1,-1};
  float ct_reset[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, ct_ignore[16]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
  FILE *fp;

   if( filename != NULL ){ len = strlen(filename);
      if( (fp=fopen(filename,"r")) == NULL ){
         fprintf(stderr,"load_config: cant open %s to read\n", filename);
         return(-1);
      }
      //fprintf(stderr,"load_config: reading %s\n", filename);
      instring = len = 0;// read line by line, copy to conf_data, skip space
      while( fgets(tmp, 80, fp) != NULL ){ // tmp always null-terminated
         for(i=0; i<strlen(tmp); i++){ // DO NOT SKIP SPACE WITHIN STRINGS
            if( tmp[i] == '"' ){ instring = 1-instring; }
            if( !isspace(tmp[i]) || instring ){ config_data[len++] = tmp[i]; }
         }
      }
      fclose(fp);
   } else if( buffer != NULL ){
      instring = len = 0; // read line by line, copy to conf_data, skip space
      for(i=0; i<strlen(buffer); i++){    // DO NOT SKIP SPACE WITHIN STRINGS
         if( buffer[i] == '"' ){ instring = 1-instring; }
         if( !isspace(buffer[i])||instring ){ config_data[len++]=buffer[i]; }
      }
   } else {
      fprintf(stderr,"load_config: no file or buffer specified\n");
   }
   clear_config(cfg); // setup hardcoded variables
   ptr=config_data;
   if( strncmp(ptr,"{\"Analyzer\":[", 13) != 0 ){
      fprintf(stderr,"load_config: err1 byte %ld\n", ptr-config_data);
       return(-1);
   } ptr += 13;
   if( strncmp(ptr,"{\"Variables\":[", 14) != 0 ){
      fprintf(stderr,"load_config: err2 byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 14;
   while( 1 ){ // variables
      if( strncmp(ptr,"]},", 3) == 0 ){ ptr+=3; break; }// empty section
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: err3 byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"title\":\"", 10) != 0 ){
         fprintf(stderr,"load_config: err4 byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 10;
      while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      // VARIABLES CURRENTLY HARDCODED
      //  - this section of config is only used for passing to viewer
      //cfg->lock=1; add_variable(cfg, name, title);  cfg->lock=0;
      ++ptr; // skip '}'
      if( *ptr++ == ',' ){ continue; }
      ptr+=2; break; // skip '},'
   }
   if( strncmp(ptr,"{\"Gates\":[", 10) != 0 ){
      fprintf(stderr,"load_config: err5 byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 10;
   while( 1 ){ // Gates
      if( strncmp(ptr,"]},", 3) == 0 ){ ptr+=3; break; }// empty section
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: err6 byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      cfg->lock=0;
      if( strncmp(ptr,",\"gateCondition\":[", 18) != 0 ){
         fprintf(stderr,"load_config: err7 byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 18;
      while( 1 ){ // Gate-Conditions
         if( strncmp(ptr,"{\"indexID\":", 11) != 0 ){
            fprintf(stderr,"load_config: err8 byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 11; while( *ptr != ',' ){ ++ptr; /* skip index id */ }
         if( strncmp(ptr,",\"Variable\":\"", 13) != 0 ){
            fprintf(stderr,"load_config: err8b byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 13; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
         if( strncmp(ptr,",\"Logic\":\"", 10) != 0 ){
            fprintf(stderr,"load_config: err8c byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 10;
         if(       strncmp(ptr,"GE",2) == 0 ){ sprintf(op,">=");
         } else if(strncmp(ptr,"GT",2) == 0 ){ sprintf(op,">");
         } else if(strncmp(ptr,"LE",2) == 0 ){ sprintf(op,"<=");
         } else if(strncmp(ptr,"LT",2) == 0 ){ sprintf(op,"<");
         } else if(strncmp(ptr,"EQ",2) == 0 ){ sprintf(op,"=");
         } else if(strncmp(ptr,"RA",2) == 0 ){ sprintf(op,"RA");
         } else {
            fprintf(stderr,"load_config:err9 byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 3;
         if( strncmp(ptr,",\"Value\":", 9) != 0 ){
            fprintf(stderr,"load_config: err9b byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 9; valstr=ptr; while( isdigit(*ptr) ){ ++ptr; } *ptr++ = 0;
         if( sscanf( valstr, "%d", &value) < 1 ){
            fprintf(stderr,"load_config:errA byte %ld\n", ptr-config_data);
            return(-1);
         }
         cfg->lock=0;
         if( *ptr++ == ',' ){ continue; }
         ++ptr; break; // skip ']'
      }
      if( *ptr++ == ',' ){ continue; }
      ptr += 2; break; // skip '},'
   }
   if( strncmp(ptr,"{\"Histograms\":[", 15) != 0 ){
      fprintf(stderr,"load_config: errB byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 15;
   while( 1 ){ // Histograms
      if( strncmp(ptr,"]},", 3) == 0 ){ ptr+=3; break; }// empty section
      if( strncmp(ptr,"{\"name\":\"",9) != 0 ){
         fprintf(stderr,"load_config: errC byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      //if( strncmp(ptr,"{\"title\":\"",10) != 0 ){
      //   fprintf(stderr,"load_config: errC byte %ld\n", ptr-config_data);
      //   return(-1);
      //} ptr += 10;
      //title = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"path\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errD byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"Xvariable\":\"", 14) != 0 ){
         fprintf(stderr,"load_config: errE byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 14;
      while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"Xmin\":", 8) != 0 ){
         fprintf(stderr,"load_config: errFa byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 8;
      valstr = ptr; while( isdigit(*ptr) ){ ++ptr; } tmp[0]=*ptr; *ptr++=0;
      if( sscanf( valstr, "%d", &val2) < 1 ){
         fprintf(stderr,"load_config:errFb byte %ld\n", ptr-config_data);
            return(-1);
      }
      if( strncmp(ptr,"\"Xmax\":", 7) != 0 ){
         fprintf(stderr,"load_config: errFc byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 7;
      valstr = ptr; while( isdigit(*ptr) ){ ++ptr; } tmp[0]=*ptr; *ptr++=0;
      if( sscanf( valstr, "%d", &val3) < 1 ){
         fprintf(stderr,"load_config:errFd byte %ld\n", ptr-config_data);
            return(-1);
      }
      if( strncmp(ptr,"\"Xbins\":", 8) != 0 ){
         fprintf(stderr,"load_config: errFe byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 8;
      valstr = ptr; while( isdigit(*ptr) ){ ++ptr; } tmp[0]=*ptr; *ptr++=0;
      if( sscanf( valstr, "%d", &value) < 1 ){
         fprintf(stderr,"load_config:errFf byte %ld\n", ptr-config_data);
            return(-1);
      }
      if( strncmp(ptr,"\"Yvariable\":\"", 13) != 0 ){
         cfg->lock=0;
      } else {
         ptr += 13;
         while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
         if( strncmp(ptr,",\"Ymin\":", 8) != 0 ){
            fprintf(stderr,"load_config: errG byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 8;
         valstr = ptr; while( isdigit(*ptr) ){++ptr;} tmp[0]=*ptr;*ptr++=0;
         if( sscanf( valstr, "%d", &val5) < 1 ){
            fprintf(stderr,"load_config:errH byte %ld\n", ptr-config_data);
               return(-1);
         }
         if( strncmp(ptr,"\"Ymax\":", 7) != 0 ){
            fprintf(stderr,"load_config: errHa byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 7;
         valstr = ptr; while( isdigit(*ptr) ){++ptr;} tmp[0]=*ptr;*ptr++=0;
         if( sscanf( valstr, "%d", &val6) < 1 ){
            fprintf(stderr,"load_config:errHb byte %ld\n", ptr-config_data);
               return(-1);
         }
         if( strncmp(ptr,"\"Ybins\":", 8) != 0 ){
            fprintf(stderr,"load_config: errHc byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 8;
         valstr = ptr; while( isdigit(*ptr) ){++ptr;} tmp[0]=*ptr;*ptr++=0;
         if( sscanf( valstr, "%d", &val4) < 1 ){
            fprintf(stderr,"load_config:errI byte %ld\n", ptr-config_data);
               return(-1);
         }
         cfg->lock=0;
      }
      if( strncmp(ptr,"\"histogramCondition\":[", 22) != 0 ){
         fprintf(stderr,"load_config: errJ byte %ld\n", ptr-config_data);
            return(-1);
      } ptr += 22;
      while(1){  // Histo gates
         if( strncmp(ptr,"]}", 2) == 0 ){ ptr+=2; break; } // empty list
         if( strncmp(ptr,"{\"indexID\":", 11) != 0 ){
            fprintf(stderr,"load_config: errK byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 11; while( isdigit(*ptr) || *ptr=='-' ){++ptr;}
         if( strncmp(ptr,",\"Gate\":\"", 9) != 0 ){
            fprintf(stderr,"load_config: errKa byte %ld\n", ptr-config_data);
            return(-1);
         } ptr += 9;
         valstr = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
         ++ptr; // skip '}' - end of single condition
         cfg->lock=0;
         if( *ptr == ',' ){ ++ptr; continue; } else { ptr +=2; break;}// ']}'
      }
      if( strncmp(ptr,",{\"", 3) == 0 ){ ++ptr; continue; }
      ptr += 3; break; // skip closing ]},
   }
   if( strncmp(ptr,"{\"Globals\":[", 12) != 0 ){
      fprintf(stderr,"load_config: errL byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 12;
   while( 1 ){ // Globals
      if( strncmp(ptr,"]},", 3) == 0 ){ ptr+=3; break; }// empty section
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errM byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"min\":",7) != 0 ){
         fprintf(stderr,"load_config: errN byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 7;
      valstr = ptr; while( isdigit(*ptr) || *ptr=='-' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%d", &value) < 1 ){
         fprintf(stderr,"load_config:errO byte %ld\n", ptr-config_data);
         return(-1);
      }
      if( strncmp(ptr,"\"max\":",6) != 0 ){
         fprintf(stderr,"load_config: errP byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 6;
      valstr = ptr; while( isdigit(*ptr) || *ptr=='-' ){++ptr;} *ptr++ = 0;
      // last char written over with 0 was '}'
      if( sscanf( valstr, "%d", &val2) < 1 ){
         fprintf(stderr,"load_config:errQ byte %ld\n", ptr-config_data);
         return(-1);
      }
      cfg->lock=0;
      if( *ptr++ == ',' ){ continue; } // have skipped ']' if not
      ptr+=2; break; // skip '},'
   }
   if( strncmp(ptr,"{\"Calibrations\":[", 17) != 0 ){
      fprintf(stderr,"load_config: errR byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 17;
   while( 1 ){ // Calibrations
      if( strncmp(ptr,"]},", 3) == 0 ){ ptr+=3; break; }// empty section or end of section
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errS byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"address\":",11) != 0 ){
         fprintf(stderr,"load_config: errV byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 11; valstr = ptr;
      while( isdigit(*ptr) || *ptr=='-' || *ptr=='.' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%d", &address) < 1 ){
         fprintf(stderr,"load_config:errX byte %ld\n", ptr-config_data);
         return(-1);
      }
      if( strncmp(ptr,"\"datatype\":",11) != 0 ){
         fprintf(stderr,"load_config: errV byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 11; valstr = ptr;
      while( isdigit(*ptr) || *ptr=='-' || *ptr=='.' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%d", &type) < 1 ){
        fprintf(stderr,"load_config:errX byte %ld\n", ptr-config_data);
        return(-1);
      }
      if( strncmp(ptr,"\"offset\":",9) != 0 ){
        fprintf(stderr,"load_config: errT byte %ld\n", ptr-config_data);
        return(-1);
      } ptr += 9; valstr = ptr;
      while( isdigit(*ptr) || *ptr=='-' || *ptr=='.' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%f", &offset) < 1 ){
        fprintf(stderr,"load_config:errU byte %ld\n", ptr-config_data);
        return(-1);
      }
      if( strncmp(ptr,"\"gain\":",7) != 0 ){
        fprintf(stderr,"load_config: errV byte %ld\n", ptr-config_data);
        return(-1);
      } ptr += 7; valstr = ptr;
      while( isdigit(*ptr) || *ptr=='-' || *ptr=='.' ){++ptr;} *ptr++ = 0;
      if( sscanf( valstr, "%f", &gain) < 1 ){
        fprintf(stderr,"load_config:errX byte %ld\n", ptr-config_data);
        return(-1);
      }
      if( strncmp(ptr,"\"quad\":",7) != 0 ){
        fprintf(stderr,"load_config: errY byte %ld\n", ptr-config_data);
        return(-1);
      } ptr += 7; valstr = ptr;
      while(isdigit(*ptr)||*ptr=='.'||*ptr=='-'||*ptr=='+'||*ptr=='e'){++ptr;}
      *ptr++ = 0; // last char written over with 0 was '}'
      if( sscanf( valstr, "%f", &quad) < 1 ){
        fprintf(stderr,"load_config:errZ byte %ld\n", ptr-config_data);
        return(-1);
      }
      // The pileup correction parameters were introduced in Feb 2025.
      // Config files before this date will not have pileup corrections, and after this they are optional
      if( strncmp(ptr,"\"pileupk1\":",11) == 0 ){
        ptr += 11; valstr = ptr;
        if( sscanf(valstr, "[%f,%f,%e,%e,%e,%e,%e], ", &puk1[0],&puk1[1],&puk1[2],&puk1[3],&puk1[4],&puk1[5],&puk1[6]) != 7 ){
          fprintf(stderr,"load_config:errPUA byte %ld\n", ptr-config_data);
          return(-1);
        }
        while(*ptr!='\"'){++ptr;}
        if( strncmp(ptr,"\"pileupk2\":",11) != 0 ){
          fprintf(stderr,"load_config:errPUB byte %ld\n", ptr-config_data);
          return(-1);
        } ptr += 11; valstr = ptr;
        if( sscanf(valstr, "[%f,%f,%e,%e,%e,%e,%e], ", &puk2[0],&puk2[1],&puk2[2],&puk2[3],&puk2[4],&puk2[5],&puk2[6]) != 7 ){
          fprintf(stderr,"load_config:errPUC byte %ld\n", ptr-config_data);
          return(-1);
        }
        while(*ptr!='\"'){++ptr;}
        if( strncmp(ptr,"\"pileupE1\":",11) != 0 ){
          fprintf(stderr,"load_config:errPUD byte %ld\n", ptr-config_data);
          return(-1);
        } ptr += 11; valstr = ptr;
        if( sscanf(valstr, "[%f,%f,%e,%e,%e,%e,%e]", &puE1[0],&puE1[1],&puE1[2],&puE1[3],&puE1[4],&puE1[5],&puE1[6]) != 7 ){
          fprintf(stderr,"load_config:errPUE byte %ld\n", ptr-config_data);
          return(-1);
        }
        while(*ptr!=']'){++ptr;}
        ptr+=2;
      }else if(strncmp(name,"GRG",3)==0){ // Only process pileup and crosstalk for HPGe
           memcpy(puk1,puk_reset, 7 * sizeof(float));
           memcpy(puk2,puk_reset, 7 * sizeof(float));
           memcpy(puE1,puE1_reset, 7 * sizeof(float));
      }else{
            memcpy(puk1,pu_ignore, 7 * sizeof(float));
            memcpy(puk2,pu_ignore, 7 * sizeof(float));
            memcpy(puE1,pu_ignore, 7 * sizeof(float));
      }
      // The crosstalk correction parameters were introduced in June 2025.
      // Config files before this date will not have crosstalk corrections, and after this they are optional
      if( strncmp(ptr,"\"crosstalk0\":",13) == 0 ){
        ptr += 13; valstr = ptr;
        if( sscanf(valstr, "[%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f], ", &ct0[0],&ct0[1],&ct0[2],&ct0[3],&ct0[4],&ct0[5],&ct0[6],&ct0[7],&ct0[8],&ct0[9],&ct0[10],&ct0[11],&ct0[12],&ct0[13],&ct0[14],&ct0[15]) != 16 ){
          fprintf(stderr,"load_config:errCTA byte %ld\n", ptr-config_data);
          return(-1);
        }
        while(*ptr!='\"'){++ptr;}
        if( strncmp(ptr,"\"crosstalk1\":",13) != 0 ){
          fprintf(stderr,"load_config:errCTB byte %ld\n", ptr-config_data);
          return(-1);
        } ptr += 13; valstr = ptr;
        if( sscanf(valstr, "[%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f], ", &ct1[0],&ct1[1],&ct1[2],&ct1[3],&ct1[4],&ct1[5],&ct1[6],&ct1[7],&ct1[8],&ct1[9],&ct1[10],&ct1[11],&ct1[12],&ct1[13],&ct1[14],&ct1[15]) != 16 ){
          fprintf(stderr,"load_config:errCTC byte %ld\n", ptr-config_data);
          return(-1);
        }
        while(*ptr!='\"'){++ptr;}
        if( strncmp(ptr,"\"crosstalk2\":",13) != 0 ){
          fprintf(stderr,"load_config:errCTD byte %ld\n", ptr-config_data);
          return(-1);
        } ptr += 13; valstr = ptr;
        if( sscanf(valstr, "[%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f]", &ct2[0],&ct2[1],&ct2[2],&ct2[3],&ct2[4],&ct2[5],&ct2[6],&ct2[7],&ct2[8],&ct2[9],&ct2[10],&ct2[11],&ct2[12],&ct2[13],&ct2[14],&ct2[15]) != 16 ){
          fprintf(stderr,"load_config:errCTE byte %ld\n", ptr-config_data);
          return(-1);
        }
        while(*ptr!=']'){++ptr;}
        ptr+=2;
      }
      cfg->lock=1; edit_calibration(cfg,name,offset,gain,quad,puk1,puk2,puE1,ct0,ct1,ct2,address,type,1); cfg->lock=0;
      if( *ptr++ == ',' ){ continue; } // have skipped ']' if not
      ptr+=2; break; // skip '},'
    }
   if( strncmp(ptr,"{\"Directories\":[", 16) != 0 ){
      fprintf(stderr,"load_config: errZA byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 16;
   while(1){
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errZB byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"Path\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errZC byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      cfg->lock=0;
      ++ptr; // skip '}'
      if( *ptr++ == ',' ){ continue; }
      ptr+=2; break; // skip '},'
   }
   if( strncmp(ptr,"{\"Midas\":[", 10) != 0 ){
      fprintf(stderr,"load_config: errZD byte %ld\n", ptr-config_data);
      return(-1);
   } ptr += 10;
   while(1){
      if( strncmp(ptr,"{\"name\":\"", 9) != 0 ){
         fprintf(stderr,"load_config: errZE byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 9;
      name = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      if( strncmp(ptr,",\"Value\":\"", 10) != 0 ){
         fprintf(stderr,"load_config: errZF byte %ld\n", ptr-config_data);
         return(-1);
      } ptr += 10;
      valstr = ptr; while( *ptr != '"' ){ ++ptr; } *ptr++ = 0;
      cfg->lock=1; set_midas_param(cfg, name, valstr);  cfg->lock=0;
      ++ptr; // skip '}'
      if( *ptr++ == ',' ){ continue; }
      ++ptr; break; // skip '}'
   }
   if( strncmp(ptr,"]}", 2) != 0 || ptr+2-config_data != len ){
      fprintf(stderr,"load_config: errR near %ld\n", ptr-config_data);
      return(-1);
   }
   return(0);
}

///////////////////////////////////////////////////////////////////////////
/////////////////          Config Commands          ///////////////////////
///////////////////////////////////////////////////////////////////////////
// config changes during sorting ... (sort probably running most of time)
// - sort will need to grab config and copy so any changes do not affect it
//


// init_config() call on server startup ...
//   get most recent config and initialise all variables, gates, histos
// *In case no recent config exists, a default config is setup first*
//  (this is usually immediately overwritten by the recent config)
//  (could change this to only setup the default, if recent fails)
int init_config(const char *filename)
{
   Config *cfg = configs[0];
   if( cfg == NULL ){ // not yet alloc'd live set
      if( (cfg=configs[0]=add_config("live")) == NULL ){ return(-1); }
      if( (configs[1]=add_config("sort")) == NULL ){ return(-1); }
      configs[0]->type = configs[1]->type = MEM_CONFIG;
   }
   load_config(cfg, filename, NULL); // attempt to load, ignore any error
   strncpy(cfg->configName,filename,SYS_PATH_LENGTH-1);
   //clear_calibrations(cfg); // Clear the calibrations to default values following server restart
   fprintf(stdout,"Configuration loaded from file: %s\n",filename);
   return(0);
}

int clear_config(Config *cfg)
{
   time_t current_time = time(NULL);
   int i;
   memset(cfg, 0, sizeof(Config));      // delete any current vars etc.
   for(i=0; i<MAX_CALIB;     i++){cfg->calib[i]    = &cfg->calib_array[i]; }
   for(i=0; i<MAX_GLOBALS;   i++){cfg->globals[i]  = &cfg->global_array[i]; }
   for(i=0; i<MAX_CONDS;     i++){cfg->condlist[i] = &cfg->cond_array[i]; }
   for(i=0; i<MAX_GATES;     i++){cfg->gatelist[i] = &cfg->gate_array[i]; }
   // used_sortvars, and user histos start empty and are filled randomly
   cfg->mtime = current_time;
   return(0);
}

int clear_calibrations(Config *cfg)
{
  float offset=0, gain=1, quad=0;
  int i, address=-1, datatype=-1;
  // Initialize values to defaults
  // Values of -1 are ignored by edit_calibration - use this for all channels that are not HPGe to avoid bloating the size of the config
  float puk_reset[7]={1,0,0,0,0,0,0}, puE1_reset[7]={0,0,0,0,0,0,0}, pu_ignore[7]={-1,-1,-1,-1,-1,-1,-1};
  float ct_reset[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, ct_ignore[16]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

  // delete any calibration values
  for(i=0; i<cfg->ncal;     i++){
    if(strncmp(cfg->calib[i]->name,"GRG",3)==0){
      edit_calibration(cfg, cfg->calib[i]->name, offset, gain, quad, puk_reset, puk_reset, puE1_reset, ct_reset, ct_reset, ct_reset, address, datatype, 1);
    }else{
      edit_calibration(cfg, cfg->calib[i]->name, offset, gain, quad, pu_ignore, pu_ignore, pu_ignore, ct_ignore, ct_ignore, ct_ignore, address, datatype, 1);
    }
  }
  fprintf(stdout,"Cleared all calibrations.\n");
  return(0);
}

int copy_config(Config *src, Config *dst)
{
   int i;

   src->lock = 1;
   memset(dst, 0, sizeof(Config));      // delete any current vars, gates etc.
   memcpy(dst, src, sizeof(Config));    // add all of above from live config
   // below is wrong - src array lists can contain holes if things were deleted
   //    use same offsets as in src arrays:
   //          offset = src->calib[i] - &src->calib_array[0];
   //   dst->calib[i] = offset + &dst->calib_array[0]
   for(i=0; i<MAX_CALIB;     i++){dst->calib[i]    = &dst->calib_array[i]; }
   for(i=0; i<MAX_GLOBALS;   i++){dst-> globals[i] = &dst->global_array[i]; }
   for(i=0; i<MAX_CONDS;     i++){dst->condlist[i] = &dst->cond_array[i]; }
   for(i=0; i<MAX_GATES;     i++){dst->gatelist[i] = &dst->gate_array[i]; }

   // some of the arrays contain pointers: cond_array has var pointers
   //                                      gate_array has cond pointers
   // these have to be copied the long way
   dst->nconds=0;
   dst->ngates=0;

   // copy config histograms  ODB histos will follow later
   /*for(i=0; i<src->nhistos; i++){
      histo = src->histo_list[i];
      tmp = ( histo->ybins ) ? histo->yvar->name : NULL;
      if( add_histo(dst, histo->handle, histo->title, histo->path, histo->xbins, histo->xvar->name, 0, histo->xbins, histo->ybins, tmp, 0, histo->ybins) ){ return(-1); }
      // apply gates ...
      for(j=0; j<histo->num_gates; j++){
         apply_gate(dst, histo->handle, histo->gate_names[j]);
      }
   }*/
   /*
   // update user_histo_list to point to new histos
   for(i=0; i<src->nuser; i++){
      srchist = src->user_histos[i];
      if( (histo = find_histo(dst, srchist->handle) ) == NULL ){
            printf("copy_config: impossible error#1\n"); continue;
      }
      dst->user_histos[i] = histo;
   }
   // update used_vars list to point to dst->sortvars (from src->sortvars)
   for(i=0; i<src->nusedvar; i++){
      srcvar = src->usedvars[i];
      if( (dstvar = find_sortvar(dst, srcvar->name)) == NULL ){
         printf("copy_config: impossible error#2\n"); continue;
      }
      dst->usedvars[i] = dstvar;
   }
   // update sortvar histo list to point to dst histos
   for(i=0; i<src->nsortvar; i++){
      srcvar = &src->varlist[i];
      dstvar = &dst->varlist[i];
      for(j=0; j<srcvar->use_count_x; j++){
         srchist = srcvar->histo_list_x[i];
         if( (histo = find_histo(dst, srchist->handle) ) == NULL ){
            printf("copy_config: impossible error#3\n"); continue;
         }
         dstvar->histo_list_x[i] = histo;
      }
   }
   */
   /* tmp  = (char *)dst;
   tmp2 = (char *)src;
   for(i=0; i<sizeof(Config); i+=8){
      tmp3 = *(long *)(tmp+i);
      ptr = (char *)tmp3;
      if( ptr >= tmp2 && ptr <= tmp2+sizeof(Config) ){
         printf("offset:%d\n", i);
      }
   }
   */
   //dst->odb_daqsize = src->odb_daqsize;
   src->lock = 0; return(0);
}

Config *add_config(char *name)
{
   Config *cfg;
   int i, len;

   for(i=0; i<MAX_CONFIGS; i++){
      if( configs[i] == NULL ){ break; }
   }
   if( i == MAX_CONFIGS ){
      fprintf(stderr,"Exceeded Max histogram sets\n");
   }
   if( (cfg = configs[i] = calloc(1, sizeof(Config))) == NULL ){
      fprintf(stderr,"can't alloc new config\n");
   }
   if( (len=strlen(name)+1) >  SYS_PATH_LENGTH ){
      fprintf(stderr,"truncating configname: %s\n", name);
      len = SYS_PATH_LENGTH;
   }
   memcpy(cfg->name, name, len);
   for(i=0; i<MAX_CALIB;   i++){cfg->calib[i]    = &cfg->calib_array[i]; }
   for(i=0; i<MAX_GLOBALS; i++){cfg->globals[i] = &cfg->global_array[i]; }
   for(i=0; i<MAX_CONDS;   i++){cfg->condlist[i] = &cfg->cond_array[i]; }
   for(i=0; i<MAX_GATES;  i++){cfg->gatelist[i] = &cfg->gate_array[i]; }
   return( cfg );
}

int remove_config(Config *cfg)
{
   int i;

   if( cfg == configs[0] || cfg == configs[0] ){ return(0); }// do not remove set#0 or #1
   for(i=1; i<MAX_CONFIGS; i++){
      if( cfg == configs[i] ){ break; }
   }
   if( i == MAX_CONFIGS ){
      fprintf(stderr,"can't find config to remove\n");
   }
   clear_config(cfg);
   free(cfg);
   configs[i] = NULL;
   return(0);
}
/////////////////////////////   CALIBRATION   /////////////////////////////

int edit_calibration(Config *cfg, char *name, float offset, float gain, float quad, float puk1[7], float puk2[7], float puE1[7], float ct0[16], float ct1[16], float ct2[16], int address, int type, int overwrite)
{
   time_t current_time = time(NULL);
   int i,j, len;
   Cal_coeff *cal;

   for(i=0; i<cfg->ncal; i++){ cal = cfg->calib[i];
      if( strncmp(name, cfg->calib[i]->name, strlen(name)) == 0 &&
          strlen(name) == strlen(cfg->calib[i]->name)       ){ break; }
   }
   if( i < cfg->ncal ){ // calib already exists
      if( overwrite ){
         if( offset != -1 ){ cal->offset = offset; }
         if( gain   != -1 ){ cal->gain = gain; }
         if( quad   != -1 ){ cal->quad = quad; }
         if( puk1[0] != -1 ){ for(j=0; j<7; j++){cal->pileupk1[j] = puk1[j];} }
         if( puk2[0] != -1 ){ for(j=0; j<7; j++){cal->pileupk2[j] = puk2[j];} }
         if( puE1[0] != -1 ){ for(j=0; j<7; j++){cal->pileupE1[j] = puE1[j];} }
         if( ct0[0] != -1 ){ for(j=0; j<16; j++){cal->crosstalk0[j] = ct0[j];} }
         if( ct1[0] != -1 ){ for(j=0; j<16; j++){cal->crosstalk1[j] = ct1[j];} }
         if( ct2[0] != -1 ){ for(j=0; j<16; j++){cal->crosstalk2[j] = ct2[j];} }
      }
      if( address != -1 ){
         cal->address = address;  cal->datatype = type;
      }
   } else { // not already there ... add new one
      if( cfg->ncal >= MAX_CALIB ){
         fprintf(stderr,"too many calibs: %d\n", MAX_CALIB); return(-1);
      }
      cal = cfg->calib[i]; ++cfg->ncal;
      if( (len=strlen(name)+1) > 64 ){
         fprintf(stderr,"truncating calibname: %s\n", name);
         len = 64;
      }
      memcpy(cal->name, name, len);
      if( offset != -1 ){ cal->offset = offset; }else{ cal->offset = 0; }
      if( gain   != -1 ){ cal->gain = gain; }else{ cal->gain = 1; }
      if( quad   != -1 ){ cal->quad = quad; }else{ cal->quad = 0; }
      if( puk1[0] != -1 ){
        for(j=0; j<7; j++){ cal->pileupk1[j] = puk1[j]; }
      }else{
        cal->pileupk1[0]=1; cal->pileupk1[1]=0; cal->pileupk1[2]=0; cal->pileupk1[3]=0; cal->pileupk1[4]=0; cal->pileupk1[5]=0; cal->pileupk1[6]=0;
      }
      if( puk2[0] != -1 ){
        for(j=0; j<7; j++){cal->pileupk2[j] = puk2[j]; }
      }else{
        cal->pileupk2[0]=1; cal->pileupk2[1]=0; cal->pileupk2[2]=0; cal->pileupk2[3]=0; cal->pileupk2[4]=0; cal->pileupk2[5]=0; cal->pileupk2[6]=0;
      }
      if( puE1[0] != -1 ){
        for(j=0; j<7; j++){cal->pileupE1[j] = puE1[j]; }
      }else{
        cal->pileupE1[0]=0; cal->pileupE1[1]=0; cal->pileupE1[2]=0; cal->pileupE1[3]=0; cal->pileupE1[4]=0; cal->pileupE1[5]=0; cal->pileupE1[6]=0;
      }
      if( ct0[0] != -1 ){
        for(j=0; j<16; j++){ cal->crosstalk0[j] = ct0[j]; }
      }else{
        cal->crosstalk0[0]=0; cal->crosstalk0[1]=0; cal->crosstalk0[2]=0; cal->crosstalk0[3]=0; cal->crosstalk0[4]=0; cal->crosstalk0[5]=0; cal->crosstalk0[6]=0; cal->crosstalk0[7]=0;
        cal->crosstalk0[8]=0; cal->crosstalk0[9]=0; cal->crosstalk0[10]=0; cal->crosstalk0[11]=0; cal->crosstalk0[12]=0; cal->crosstalk0[13]=0; cal->crosstalk0[14]=0; cal->crosstalk0[15]=0;
      }
      if( ct1[0] != -1 ){
        for(j=0; j<16; j++){ cal->crosstalk1[j] = ct1[j]; }
      }else{
        cal->crosstalk1[0]=0; cal->crosstalk1[1]=0; cal->crosstalk1[2]=0; cal->crosstalk1[3]=0; cal->crosstalk1[4]=0; cal->crosstalk1[5]=0; cal->crosstalk1[6]=0; cal->crosstalk1[7]=0;
        cal->crosstalk1[8]=0; cal->crosstalk1[9]=0; cal->crosstalk1[10]=0; cal->crosstalk1[11]=0; cal->crosstalk1[12]=0; cal->crosstalk1[13]=0; cal->crosstalk1[14]=0; cal->crosstalk1[15]=0;
      }
      if( ct2[0] != -1 ){
        for(j=0; j<16; j++){ cal->crosstalk2[j] = ct2[j]; }
      }else{
        cal->crosstalk2[0]=0; cal->crosstalk2[1]=0; cal->crosstalk2[2]=0; cal->crosstalk2[3]=0; cal->crosstalk2[4]=0; cal->crosstalk2[5]=0; cal->crosstalk2[6]=0; cal->crosstalk2[7]=0;
        cal->crosstalk2[8]=0; cal->crosstalk2[9]=0; cal->crosstalk2[10]=0; cal->crosstalk2[11]=0; cal->crosstalk2[12]=0; cal->crosstalk2[13]=0; cal->crosstalk2[14]=0; cal->crosstalk2[15]=0;
      }
      cal->address = address;  cal->datatype = type;
    }
    //printf("saving config edit_calibration\n");
   cfg->mtime = current_time;
   return(0);
}

int set_directory(Config *cfg, char *name, char *path)
{
   time_t current_time = time(NULL);
   int len;


   if( (len=strlen(path)) >= SYS_PATH_LENGTH ){
      fprintf(stderr,"set_directory: path too long[%s]\n", path);
      return(-1);
   }
   if( strncmp(name, "Data",   4) == 0 ){
      memcpy(cfg->data_dir, path, len+1);
   } else if( strncmp(name, "Config", 6) == 0 ){
      memcpy(cfg->config_dir, path, len+1);
   } else {
      fprintf(stderr,"set_directory: Unknown directory:%s\n", name);
      return(-1);
   }
   //printf("saving config set_directory\n");
   cfg->mtime = current_time;
   return(0);
}

int set_midas_param(Config *cfg, char *name, char *value)
{
   time_t current_time = time(NULL);
   int len;
   char clean_string[128], *tmp;

   if( (len=strlen(value)) >= SYS_PATH_LENGTH ){
      fprintf(stderr,"set_midas_param: value too long[%s]\n", value);
      return(-1);
   }
   if(strncmp(name,"Title", 5) == 0 ){
     sprintf(clean_string, "%s", value);
     if( (tmp = strstr(clean_string, "\t"))>0 ){ // This illegal character cannot be handled in the browser
         while( (tmp = strstr(clean_string, "\t"))>0 ){
            strncpy(tmp, " ", 1); // keep the length the same, and make use of the terminating character already in clean_string
         }
         memcpy(cfg->midas_title, clean_string, len+1);
      }else{
         memcpy(cfg->midas_title, value, len+1);
      }
   } else if( strncmp(name, "StartTime",  9) == 0 ){
      if( sscanf(value, "%d", &cfg->midas_start_time) < 1 ){
         fprintf(stderr,"set_midas_param: can't read starttime: %s\n", value);
      }
   } else if( strncmp(name, "Duration", 8) == 0 ){
      if( sscanf(value, "%d", &cfg->midas_runtime) < 1 ){
         fprintf(stderr,"set_midas_param: can't read runtime: %s\n", value);
      }
   } else {
      fprintf(stderr,"set_midas_param: Unknown param:%s\n", name);
      return(-1);
   }
   //printf("saving config set_midas_param\n");
   cfg->mtime = current_time;
   return(0);
}


#define READ_LIMIT (5*1024*1024) // don't read more than this looking for odb
int read_datafile_info(Sort_status *sort, char *path)
{
   int i, expt=0, ppg=0, flt=0, done=0;
   char tmp[256], *ptr;
   FILE *fp;
   if( (fp=fopen(path,"r")) == NULL ){
      fprintf(stderr,"read_datafile_info: can't read file: %s\n", path);
      return(-1);
   }
   while( (ptr=fgets(tmp, 256, fp)) != NULL ){
      if( ftell(fp) > READ_LIMIT ){ break; }
      while( isspace(*ptr) ){ ++ptr; }
      if( strncmp(ptr,"</dir>",6) == 0 ){ expt = 0; }
      if( strncmp(ptr,"<dir name=\"Run parameters\">",27) == 0 ){ expt = 2; }
      // The Title/comment are in Run parameters
      //    Experiment/Status-items and Edit-on-start contain LINKS
      //if( strncmp(ptr,"<dir name=\"Experiment\">",23)     == 0 ){ expt = 1; }
      //if( strncmp(ptr,"<dir name=\"Edit on start\">", 26) == 0 ){ expt = 3; }
      if( expt && strncmp(ptr,"<key name=\"Run Title\"",21) == 0 ){
         ptr+=21; while(*ptr != '>'){ ++ptr; } ++ptr; i=0;
         while(*ptr != '<'){ sort->file_info[0][i++] = *ptr++; }
         sort->file_info[0][i]=0; done |= 1;
      }
      if( expt && strncmp(ptr,"<key name=\"Comment\"",19) == 0 ){
         ptr+=19; while(*ptr != '>'){ ++ptr; } ++ptr; i=0;
         while(*ptr != '<'){ sort->file_info[1][i++] = *ptr++; }
         sort->file_info[1][i]=0; done |= 2;
      }
      if( strncmp(ptr,"<dir name=\"PPG\">",16)            == 0 ){ ppg = 1; }
      if( ppg && strncmp(ptr,"<key name=\"Current\"", 19) == 0 ){
         ptr+=19; while(*ptr != '>'){ ++ptr; } ++ptr; i=0;
         while(*ptr != '<'){ sort->file_info[2][i++] = *ptr++; }
         sort->file_info[2][i]=0; ppg = 0; done |= 4;
      }
      if( strncmp(ptr,"<dir name=\"Filter\">",19)         == 0 ){ flt = 1; }
      if( flt && strncmp(ptr,"<key name=\"Current\"", 19) == 0 ){
         ptr+=19; while(*ptr != '>'){ ++ptr; } ++ptr; i=0;
         while(*ptr != '<'){ sort->file_info[3][i++] = *ptr++; }
         sort->file_info[3][i]=0; flt = 0; done |= 8;
      }
      if( done == 15 ){ break; } // 1st 4 bits of done all set
   }
   fclose(fp);
   return(0);
}

int send_file_details(char *path, int fd)
{
   char tmp2[256];
   Sort_status *tmp;
   if( (tmp = calloc(sizeof(Sort_status), 1)) == NULL){
      
      fprintf(stderr,"send_file_details: failed alloc\n");
   }
   
   read_datafile_info(tmp, path);
   sprintf(tmp2,"%s\n%s\n%s\n%s\n", tmp->file_info[0], tmp->file_info[1],
                                   tmp->file_info[2], tmp->file_info[3] );
   free(tmp);
   return(0);
}

static struct stat statbuf;

// split path into name,dir  then get histo and config names
// also do stat to get size
int run_number(Sort_status *arg, char *name);
char *subrun_filename(Sort_status *sort, int subrun);
int add_sortfile(char *path)
{
   int i, plen, dlen;
   char *fname;
   Sort_status *sort;
   Config *cfg = configs[0];

   sort = get_sort_status();
   plen=strlen(path);
   for(i=plen; i>=0; i--){ if( path[i] == '/' ){ ++i; break; } }
   if( (dlen = i) == -1 ){ dlen = 0; } // no directory separator in path
   if( (sort->data_dir = malloc((size_t)((unsigned int)dlen + 2))) == NULL ){
      fprintf(stderr,"can't alloc string for data_dir");
      free_sortfile(sort); return(-1);
   }
   if( dlen == 0 ){ sprintf(sort->data_dir, ".");
   } else {
      memcpy((char *)sort->data_dir, path, dlen-1);
      *(sort->data_dir+dlen-1) = 0; // overwrite trailing '/'
      set_directory(cfg, "Data", sort->data_dir);
   }
   if( (sort->data_name = malloc(plen-dlen+1)) == NULL ){
      fprintf(stderr,"can't alloc string for :%s", path);
      free_sortfile(sort); return(-1);
   }
   memcpy((char *)sort->data_name, path+dlen, plen-dlen);
   *(sort->data_name+plen-dlen) = 0;
   if( run_number(sort, sort->data_name) ){ return(-1); }
   memset(&statbuf, 0, sizeof(struct stat)); sort->data_size = 0;

   //get number of subruns
   sort->subrun = -1;
   for(i=0; ; i++){
      fname = subrun_filename(sort, i);
         if( stat(fname, &statbuf) != 0 ){
            //fprintf(stderr,"can't stat multi-subrun: %s\n", path);
            //fprintf(stderr,"fname: %s\n", fname);
            break;
         }
         sort->data_size += (long)statbuf.st_size;
   }
   sort->num_subruns = i; sort->subrun = 0;

   return(0);
}

//////////////////////  sorting loop /////////////////////
// if( currentfile == finalfile ){ sleep(1); continue; }
// open_sortfilelist();
// sortnextfile()
//     create midas thread
//     read gains + init histos
//     sort data
//     close thread
//     write histos
// close_sortfilelist();
//////////////////////////////////////////////////////////
int open_next_sortfiles(Sort_status *sort)
{
   char tmp[256];
   if( sort->num_subruns == 0 ){
      sprintf(tmp, "%s", sort->data_name);
   } else {
      sprintf(tmp, "%s", subrun_filename(sort, 0) );
   }
   if( (sort->data_fp=fopen(tmp,"r")) == NULL ){
      fprintf(stderr,"open_next_sortfiles: can't open %s to read\n", tmp);  return(-1);
   }
   //fprintf(stdout,"sorting file: %s\n", tmp);
   sort->midas_bytes = 0;
   if( strcmp(sort->cal_src, "midas") == 0 ){
      fprintf(stdout,"Calibration method is ODB from midas file\n");
      sort->cal_overwrite = 1;
   } else {
      if( strcmp(sort->cal_src, "config") == 0 ){ 
         fprintf(stdout,"Calibration method is current config\n"); 
      }
      sort->cal_overwrite = 0;  // cal src == "config" or "file"
   }
   return(0);
}

// Only first subrun contains odb event
// *could* just sort single subrun if subrun of specified data file is nonzero
//         otherwise sort all subruns
int open_next_subrun(Sort_status *sort)
{
   char *filename;
   fclose(sort->data_fp); sort->data_fp = NULL;

   if( sort->subrun >= sort->num_subruns-1 ){ // last or no subruns
      fprintf(stderr,"Final Subrun[%d] completed\n", sort->subrun); return(-1);
   }
   filename = subrun_filename(sort, ++sort->subrun);
   if( (sort->data_fp=fopen(filename,"r")) == NULL ){
      fprintf(stderr,"open_next_subrun: can't open %s to read\n", filename);  return(-1);
   }
   fprintf(stdout,"\n================================================\n");
   fprintf(stdout,"=============  SORTING SUBRUN %3d  =============\n", sort->subrun);
   fprintf(stdout,"================================================\n\n");
   return(0);
}

char *subrun_filename(Sort_status *sort, int subrun)
{
   static char name[256];
   int len, digits;
   char tmp[64];

   //printf("run: %i, run digits: %i, subrun digits: %i\n",sort->run, sort->run_digits, sort->subrun_digits);

   sprintf(name, "%s/run", sort->data_dir);

   sprintf(tmp,"%d", sort->run); digits = strlen(tmp);
   len = strlen(name);
   while( digits++ < sort->run_digits ){ name[len] = '0'; name[++len] = 0; }
   sprintf(name+strlen(name),"%d_", sort->run);

   sprintf(tmp,"%d", subrun);  digits = strlen(tmp);
   len = strlen(name);
   while( digits++ < sort->subrun_digits ){
      name[len++]='0'; name[len] = 0;
   }
   sprintf(name+strlen(name),"%d.mid", subrun);

   return(name);
}

int close_sortfiles(Sort_status *sort)
{  // data_fp usually closed in nxtsubrun
   if( sort->online_mode ){ return(0); } // no files in this mode
   if( sort->data_fp != NULL ){ fclose(sort->data_fp); }\
   free_sortfile(sort);
   return(0);
}

int free_sortfile(Sort_status *sort)
{
   if( sort->data_dir   != NULL ){ free(sort->data_dir);   }
   if( sort->data_name  != NULL ){ free(sort->data_name);  }
   memset(sort->data_dir, 0, sizeof(sort->data_dir));
   memset(sort->data_name, 0, sizeof(sort->data_name));
   return(0);
}

// read sun/subrun from filename, then count #digits in each
int run_number(Sort_status *arg, char *name)
{
   char *ptr = name, fmt[16], tmp[256];
   FILE *fp;
   int i;

   if( strncmp(ptr, "run", 3) != 0 ){
      fprintf(stderr,"datafilename:%s does not being with \"run\"\n", name);
      return(-1);
   } ptr += 3;
   while( 1 ){
      if( !isdigit(*ptr) ){ arg->run_digits = ptr-name-3; break; }
      ++ptr;
   }
   if( *ptr != 0 ){               // name contains stuff after run number ...
      if( *ptr == '_' ){
         if( sscanf(name, "run%d_%d.mid", &arg->run, &arg->subrun) != 2 ){
            fprintf(stderr,"can't read run and subrun number in %s\n", name);
            return(-1);
         }
         arg->subrun_digits = -1; ++ptr; while( 1 ){
            if( !isdigit(*ptr) ){
               arg->subrun_digits = ptr-name-4-arg->run_digits;  break;
            }
            ++ptr;
         }
         if( arg->subrun_digits == -1 ){
            fprintf(stderr,"can't read subrun number in %s\n", name);
            return(-1);
         }
         if( strncmp(ptr, ".mid", 4) == 0 ){ return(0); }
         else {
            fprintf(stderr,"no .mid extension in datafile: %s\n", name);
            return(-1);
         }
      } else if( strncmp(ptr, ".mid", 4) == 0 ){
         fprintf(stderr,"subrun number missing in datafile: %s\n", name);
         return(-1);
      } else {
         fprintf(stderr,"bad data filename format in %s\n", name);
         return(-1);
      }
   } else {               // name only contains run number - sort all subruns
      arg->subrun = -1;
      if( sscanf(name, "run%d.mid", &arg->run) != 1 ){
         fprintf(stderr,"can't read run number in %s\n", name);
         return(-1);
      }
      for(i=1; i<5; i++){ // look for up to 5 subrun digits (usually 3)
         sprintf(fmt, "%%s/%%s_%%0%dd.mid", i);
         sprintf(tmp, fmt, arg->data_dir, name, 0 );
         if( (fp = fopen(tmp,"r")) == NULL ){ continue; }
         arg->subrun_digits = i; fclose(fp); return(0);
      }
      fprintf(stderr,"can't open subrun0 for datafile:%s\n", name);
      return(-1);
   }
   return(0);
}

int end_current_sortfile(int fd)
{
   Sort_status *arg;

   arg = get_sort_status();
   arg->end_of_data = 1; //arg-> shutdown_midas = 1;
   return(0);
}


///////////////////////////////////////////////////////////////////////////
/////////////////          Directory reading          /////////////////////
///////////////////////////////////////////////////////////////////////////
#include <dirent.h>

int send_datafile_list(char *path, int fd, int type)
{
   char tmp[256];  Sort_status *tmp_srt;
   int nlen, run, subrun, entry=0;
   struct dirent *d_ent;
   DIR *d;

   if( (d=opendir(path)) == NULL ){
      sprintf(tmp,"can't open directory %s\n",path);
      fprintf(stderr,"can't open directory %s\n", path);
      return(-1);
   }
   set_directory(configs[0], "Data", path);
   if( type == 1 ){
      if( (tmp_srt  = calloc(sizeof(Sort_status),    1)) == NULL ){
         fprintf(stderr,"send_datafile_list: failed alloc\n");
      }
   } else { tmp_srt = NULL; }
   
   while( (d_ent = readdir(d)) != NULL ){
      //fprintf(stdout,"File[%s] ...\n", d_ent->d_name);
      if( strncmp(d_ent->d_name, ".", 1) == 0 ){
         continue; // Ignore
      }
      if( sscanf(d_ent->d_name, "run%d_%d.mid", &run, &subrun) != 2 ){
         if( sscanf(d_ent->d_name, "run%d.mid", &run) != 1 ){
            continue; // Not Midas Data File
         }
         subrun=0;
      }
      nlen = strlen(d_ent->d_name);
      if( strncmp(d_ent->d_name+nlen-4, ".mid", 4)     != 0 ){ // &&
          //strncmp(d_ent->d_name+nlen-7, ".mid.gz", 7)  != 0 &&
          //strncmp(d_ent->d_name+nlen-8, ".mid.bz2", 8) != 0 ){
         continue; // Not Midas DataFilename Extension
      }
      snprintf(tmp,256,"%s/%s", path, d_ent->d_name);
      if( stat(tmp, &statbuf) != 0 ){
         fprintf(stderr,"can't stat %s\n", tmp); statbuf.st_size = 1;
      }
      //put_line(fd, d_ent->d_name, strlen(d_ent->d_name) );
      snprintf(tmp,256," , %ld ", (long)statbuf.st_size);
      //put_line(fd, tmp, strlen(tmp) );
      if( (entry % 1000) == 0 ){ printf("Entry: %d\n", entry); }
      if( type == 0 ){ continue; }

      if( subrun == 0 ){
         snprintf(tmp,256,"%s/%s", path, d_ent->d_name);
         read_datafile_info(tmp_srt, tmp);
      } else { tmp_srt->file_info[0][0] = tmp_srt->file_info[1][0] = 0; }
      if( strlen(tmp_srt->file_info[0]) > 0 ){
         snprintf(tmp,256," , %s ", tmp_srt->file_info[0] );
      } else {
         snprintf(tmp,256," , %s ", tmp_srt->file_info[1] );
      }
      //put_line(fd, tmp, strlen(tmp) );
   }
   //put_line(fd, " ]\n", 3 );
   if( tmp_srt != NULL ){ free(tmp_srt); }
   return(0);
}

