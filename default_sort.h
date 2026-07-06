// default_sort.h
// Header file for default_sort.c


//#######################################################################
//########         Subsystem and Detector definitions          ##########
//#######################################################################

// do not alter order without also changing subsys_e_vs_e, subsys_dt
#define MAX_SUBSYS      24
#define SUBSYS_HPGE_A    0
#define SUBSYS_PACES     1
#define SUBSYS_LABR_L    2
#define SUBSYS_RCMP      3
#define SUBSYS_ARIES_A   4
#define SUBSYS_ZDS_A     5 // GRIF16
#define SUBSYS_TAC_LABR  6
#define SUBSYS_LABR_BGO  7
#define SUBSYS_BGO       8
#define SUBSYS_SCEPTAR   9
#define SUBSYS_DESCANT  10
#define SUBSYS_DESWALL  11
#define SUBSYS_DSG      12
#define SUBSYS_IGNORE   13
#define SUBSYS_HPGE_B   16
#define SUBSYS_ARIES_B  17 // CAEN
#define SUBSYS_ZDS_B    18 // CAEN
#define SUBSYS_TAC_ZDS  19
#define SUBSYS_TAC_ART  20
#define SUBSYS_UNKNOWN  23
// #####################################################################

//#######################################################################
//########                PRESORT Time Gates                   ##########
//#######################################################################

// The definition of the time difference gate in 10 nanosecond units.
// The value is the maximum time difference in 10 nanosecond units.
// The default values set here are replaced by the Global value at start of sorting.

//BGO window
//300 ns (ie. 30 samples) is the GRSISort default (see https://github.com/GRIFFINCollaboration/GRSISort/blob/baf84f5947ec6a80035b01d38696b6e5d1ae2dcc/include/TAnalysisOptions.h#L70)
static int bgo_window_min = 0;
//static int bgo_window_max = 30;
static int bgo_window_max = 100;

// HPGe pileup
#define N_PU_CLASSES 15
// Pileup Class definitions
#define PU_ERROR             0
#define PU_SINGLE_HIT        1
#define PU_SINGLE_HIT_ERROR  2
#define PU_2HIT_A1ST         3
#define PU_2HIT_A2ND         4
#define PU_2HIT_B1ST         5
#define PU_2HIT_B2ND         6
#define PU_2HIT_C1ST         7
#define PU_2HIT_C2ND         8
#define PU_2HIT_ERROR        9
#define PU_3HIT_1ST         10
#define PU_3HIT_2ND         11
#define PU_3HIT_3RD         12
#define PU_3HIT_ERROR       13
#define PU_OTHER            14
