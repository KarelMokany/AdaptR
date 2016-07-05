//*****************************************************************************
// TRAFT . CPP:       (MURU)
// written by Tom Harwood Jan-Aug 2011
// Data structure to store data only for LAND, ( NULL for sea)
// and deal with multiple time steps and interpolation
// Designed to be nice and self contained so it has a realy simple interface
// Use of increment to interpolate between known data steps
// Can take allocated species data
// Also deals with dispersal or sampling over the grid
//
//*****************************************************************************
#ifndef RaftH
#define RaftH
#include "GridCell.h"
//class TGenerator;
using std::string;
 //*****************************************************************************
// TRaftGrid
// base class for data which is stored for land..
//
//*****************************************************************************
class TRaft
  {
  public:
   // constructor
   TRaft(int arg_n_rows,
		 int arg_n_cols,
		 short** arg_s_condition);
   // destructor
   ~TRaft(void);

  //GET
  /*
  long int GetNRel(void) {return n_rel;}
  long int GetNRows(void) {return n_rows;}
  long int GetNCols(void) {return n_cols;}
  int GetNLayers(void) {return n_layers;}
  int GetFromIndex(void){return i_from_stack;}
  int GetFromIndex(const long int arg_i_row,
                   const long int arg_i_col);
  int GetToIndex(void){return i_to_stack;}
  */
  //float GetGDMIntercept(void){return Run_Parameters->f_gdm_intercept;}
  TGridCell_ptr GetCell(const int arg_i_row,
                             const int arg_i_col) {if(Map[arg_i_row][arg_i_col]==NULL)return NULL; else return Map[arg_i_row][arg_i_col];}
  /*
  TGridCell_ptr GetRelCell(const int arg_i_rel,
                            const long int arg_i_row,
                            const long int arg_i_col);
  int GetRowRelCell(const int arg_i_rel,             // KM Oct 2012
                    const long int arg_i_row,
                    const long int arg_i_col);
  int GetColRelCell(const int arg_i_rel,             // KM Oct 2012
                    const long int arg_i_row,
                    const long int arg_i_col);
  float GetRelDistance(const int arg_i_rel){return f_rel_dist[arg_i_rel];}
  float GetRelBearing(const int arg_i_rel){return f_rel_bearing[arg_i_rel];}
  double GetRelProb(const int arg_i_rel){return d_rel_prob[arg_i_rel];}

  short* GetSpeciesList(const long int arg_i_row,
                        const long int arg_i_col);
  short* GetRelSpeciesList( const int arg_i_rel,
                            const long int arg_i_row,
                            const long int arg_i_col);
  //Set
  void SetDice(TGenerator* arg_Dice){Dice=arg_Dice;}
  void SetCompositionFileStream(ifstream* arg_FileStream_ptr ){CompositionFileStream=arg_FileStream_ptr;}
  int MakeLand(string arg_s_condition_file);  // Full path to load from
  int IsLand(const long int arg_i_row,
             const long int arg_i_col) {if(Map[arg_i_row][arg_i_col]==NULL)return 0; else return 1;}
  // neighbourhood stuff
  void CalculateRelativeRadialNeighbourhood(const float arg_f_radius);
  void LoadRelativeNeighbourhood(string arg_s_sample_file);
  float CalculateBearingToRel(const long int arg_i_rel_row,
                              const long int arg_i_rel_col);
  double CalculateDispersal(const float arg_f_dist);

  //Loading and managing stacks in time
  void AdvanceTime(void);
  void FlipStacks(void); // makes i_from_stack=i_to_stack and vice versa
  void InitialiseAllSeries(void);
  void RefreshAllSeries(void);

  int LoadWholeGDMStack(const int arg_i_which,  //K_FROM_STACK 0 and K_TO_STACK 1
                        string arg_s_file);
  int ReadInWholeGDMStack(const int arg_i_which,  //K_FROM_STACK 0 and K_TO_STACK 1
                          string arg_s_file);
  void InitialiseGDMSeries(void);
  void RefreshGDMSeries(void);
  int AllocateGDMStacks(void);
  void FreeGDMStacks(void);
  void MeltAllGDMStacks(void);
  int AllocateGDMIncrement(void);
  void FreeGDMIncrement(void);
  int CalculateGDMIncrement(void);

  int LoadConditionGrid(const int arg_i_which,         // K_FROM_CON 0   and  K_TO_CON  1
                        string arg_s_condition_file);  // Full path to load from
  void InitialiseCondSeries(void);
  void RefreshCondSeries(void);
  int CalculateCondIncrement(void);

  int LoadRichnessGrid(const int arg_i_which,         // K_FROM_CON 0   and  K_TO_CON  1
                        string arg_s_richness_file);  // Full path to load from
  void InitialiseRichSeries(void);
  void RefreshRichSeries(void);
  int CalculateRichIncrement(void);
  // Compositional data stuff
  void InitialiseRealisedCompositionBuffer(ifstream* arg_composition_file,
                                           const long int arg_i_start_row);
  void DestroyRealisedCompositionBuffer(void);
  void MoveRealisedCompositionBufferToNextRow(void);
  void LoadSpeciesRow(const long int arg_i_row);
  void LoadKnownSpeciesDataRow(const long int arg_i_row);
  void LoadEmptySpeciesDataRow(const long int arg_i_row);
  void FreeSpeciesRow(const long int arg_i_row);
  void FreeKnownSpeciesDataRow(const long int arg_i_row);
  */
  void LoadEmptyRaftData(int arg_n_rows,
						 int arg_n_cols,
						 int arg_n_env_variables);
  void FreeRaftData(int arg_n_rows,
					int arg_n_cols);

  int ReadInWholeENVStack(string arg_s_file,
						  int arg_n_env_variables);

  protected:
  TGridCell_ptr** Map;    // pointer to 2D array of grid cell pointers.NULL for sea, valid (pointing to a TGridCell) for land
  long int n_rows;
  long int n_cols;
  int n_layers;
  long int n_rel; // number of cells in neighbourhood ( if centre within Radius)
  short int* s_rel_row; //[17000];  // relative cell refs of neighbourhood. allocate according to requirements
  short int* s_rel_col; //[17000];  // relative cell refs of neighbourhood.
  float* f_rel_dist;
  float* f_rel_bearing;
  double* d_rel_prob;
  int i_step_years;  // number of years to subdivide the interval NB: IF ZERO then no step
  int i_step_count;  // 0 to i_step_years counter to see which year we are on COUNTS DOWN
  int i_from_stack, i_to_stack; // rotating stack references so we can juggle the from and to stacks around
  //TGenerator* Dice;
  TGridCell_ptr* AreaCells; // list of pointers to cells for the area
  int i_data_step; // flag to detect times to recalculate/ swap
  //MuruRunParameters* Run_Parameters;   // defined in structures.h
  long int* n_land_in_row; // array n_rows deep how many valid cells in each row

  // Stuff for the rolling species list
  ifstream* CompositionFileStream;    // file stream
  long int i_current_row;
  long int i_buffer_diameter;
  long int i_buffer_radius;
  long int i_min_row;
  long int i_max_row;
  }; // end class TRaft

//---------------------------------------------------------------------------
#endif
