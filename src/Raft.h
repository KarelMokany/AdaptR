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

  
  TGridCell_ptr GetCell(const int arg_i_row,
                             const int arg_i_col) {if(Map[arg_i_row][arg_i_col]==NULL)return NULL; else return Map[arg_i_row][arg_i_col];}

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
