/****************************************************************
 *
 * DEFINES.H   Written May 2011 by Tom Harwood.
 *
 * Contains constants
 *
 **********************************************************/
#ifndef TDefinesH
#define TDefinesH

//UNIX OR WINDOWS?  (\\ or /
//#define K_FOLDER "\\"  //WINDOWS
#define K_FOLDER "/"  // UNIX
//GRID DEFINITIONS
#define K_I_NO_DATA -9999
#define K_NO_DATA_ISH -500
#define K_SHORT_SCALE 5000   // Critical to relationship with compacted files
#define K_TRANS_MAX 13.125    // actually 13.127 but better safe than sorry!. Max value we can get in a unsigned short
#define K_RICH_SCALE 100   // Critical to relationship with compacted files
#define K_RICH_MAX 327    // Max value we can get in a unsigned short
#define K_SHORT_MAX 32700 // how big is a short?#define K_SHORT_MAX 32700 // how big is a short?

// Raft indices

#define K_FROM_STACK 0
#define K_TO_STACK 1
#define K_FREEZE 2
#define K_MELT -9999

// si_data_ptr[3] stack indices
#define K_NUM_DATA 7  // number of data arrays below, 0, 1 condition   2,3 relative pairwise indices
#define K_FROM_CON 0 // the from condition grid
#define K_TO_CON 1   // the to condition grid
#define K_N_S 2      // pairwise neighbour similarity  above
#define K_E_W 3      // pairwise neighbour similarity  to side
#define K_BVALUE 4   // the BIODIVERSITY stuff
#define K_NEXT_BVALUE 5  // the BIODIVERSITY stuff for the next step
#define K_ECO 6      // the eco category: Realm. Bioregion Ecoregion

#define K_MAX_BVALUE 30000 // max value of biodiversity in a cell

#define K_FROM_SURFACE 0   // calculate surface as FROM FROM
#define K_TO_SURFACE 1     // calculate surface as TO TO

#define K_UNIFORM 0
#define K_RANDOM 1

//File stuff
#define K_FLOAT 1
#define K_USHORT 0
#define K_EOF  -32000
//GENERAL STUFF
#define TRUE 1
#define FALSE 0
#define OK 1
#define BAD 0
#define K_NO_RESULT  99999999999999
#define K_NODATA -9999
#define ZERO    0
#define ONE     1
#define TWO     2
#define THREE   3
#define FOUR    4
#define SIX     6
#define TWELVE 12
#define MINUS_ONE -1

#define NOT_V_MUCH  0.001
#define NOT_V_MINUS -0.001
#define NEARLY_ZERO  0.0000000001
#define K_REALLY_BIG 99999999999999
// PI
#define PI      3.14159265358978993897
#define TWO_PI  6.28318530717957959750

// Karels stuff
#define TEN 10
#define K_DIAGONAL_WEIGHT 0.7
#define K_TWENTY_KMS 20000
#define K_TEN_KMS 10000
#define K_ONE_KM 1000
#define K_TWO_KM 2000
#define K_FOUR_KM 4000

//Conversion
#define K_HA_PER_KMSQ 100   //number of ha in a square km

#endif

