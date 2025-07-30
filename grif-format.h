#define MAX_SAMPLE_LEN   4096
#define MAX_SCALAR_LEN   256
#define PTR_BUFSIZE      4096 //changed from grif-replay
#define MAX_PSD_VALS     16

#include <stdint.h> //allows uint8_t and similiar types

// **************************************************************************
// MANY VALUES IN THIS STRUCTURE HAVE TO BE ACCESSED USING A HARDCODED OFFSET
// FROM THE START OF THE STRUCTURE - DO NOT CHANGE THE ORDERING OR
// INSERT NEW VALUES, WITHOUT ALSO ADJUSTING THE OFFSETS IN USER_SORT.C
// **************************************************************************
typedef struct griffin_fragment_struct { // was 74 bytes, now ?            //OFFSET
   long        ts;   uint16_t   address; short  deadtime;                  //0
   char      dtype;  char    array_posn; char       nhit;   char  pileup;  //3
   int          q1;  int         integ1; int          q2;   int   integ2;  //4
   int          q3;  int         integ3; int          q4;   int   integ4;  //8
   int         cfd;  int       trig_req; int    trig_acc;   int   net_id;  //12
   int   master_id;  int master_pattern; int         psd;   int cc_short;  //16
// items below are derived from items above ...
   float      ecal;  int           chan; int      subsys;   int suppress;  //20
   float      esum;  int   multiplicity; int     delta_t;   int alt_chan;  //24
   float  alt_ecal;  int            tof; int    pu_class;                  //28
} Grif_event;

// Q* are the original Q-Sum from the digitizer
// Energy are the Q-Sums divided by integ-times
// Ecal is calibrated energy#1

// ppg_pattern is now wave_expected, num_pileup

// midas-timestamp should be redundant and equal to BOR time+timestamp
//   can just check this in midas part

extern uint64_t process_event(Grif_event *ptr, int slot, FILE *out);
//extern uint64_t build_event(Grif_event *ptr, int slot, FILE *out);
extern int apply_gains(Grif_event *ptr);
extern uint64_t insert_presort_win(Grif_event *ptr, int slot, FILE *out);
extern uint64_t insert_sort_win(Grif_event *ptr, int slot, FILE *out);
extern int GetIDfromAddress(unsigned short addr);
extern int pre_sort_enter(int start_idx, int frag_idx);
extern int pre_sort_exit(int frag_idx, int end_idx);

