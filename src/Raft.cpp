//*****************************************************************************
// TRAFT . CPP:       (MURU)
// written by Tom Harwood Jan-Aug 2011
// Data structure to store data only for LAND, ( NULL for sea)
// and deal with multiple time steps and interpolation
// Designed to be nice and self contained so it has a really simple interface
// Use of increment to interpolate between known data steps
// Can take allocated species data
// Also deals with dispersal or sampling over the grid
//
//*****************************************************************************
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include "defines.h"
#include "DynamicArray.h"
#include "GridCell.h"  // this contains a lot of the stuff we use here
#include "Raft.h"
using namespace std;
//****************************************************************************
// T RAFT: Constructor
// Creates an empty Raft, which is a grid nrows*ncols of TGridCell_ptr, all
// of which are set to NULL for sea/nodata for now
// We can then upload stuff to bits of it as required
//****************************************************************************
TRaft::TRaft(int arg_n_rows,
			 int arg_n_cols,
			 short** arg_s_condition)

{
  long int i_row; // row counter
  long int i_col;// col counter
  int i_land_count;
  //string s_land_file;
  //string s_sample_file;
  //string s_file;
  n_rows=arg_n_rows;
  n_cols=arg_n_cols;

  //Run_Parameters=arg_Run_Parameters;
  //n_rows=Run_Parameters->n_rows;
  //n_cols= Run_Parameters->n_cols;
  //n_layers=Run_Parameters->n_layers;
  // Step 1: Allocate the 2D array of TGridCellpointers
  // Allocate memory for array of elements of column ..
  // i.e. list of pointers to first element of each row
  Map = new TGridCell_ptr*[arg_n_rows];
  // Allocate memory for array of elements of all the rows
  // Do this as one big block to speed things up
  TGridCell_ptr* TCD_row_ptr = new TGridCell_ptr[arg_n_rows*arg_n_cols];
  // Now point the pointers in the right place
  for(i_row=0; i_row<arg_n_rows; i_row++)
    {
	*(Map+i_row)=TCD_row_ptr; // Could also be written as T_new_pptr[i_row]=
	TCD_row_ptr+=arg_n_cols;       // skip forward to the beginning of the next row
    }
  // Allocated.. Step 2.. Nullify
  for(i_row=0;i_row<arg_n_rows;i_row++)
	{
	for(i_col=0;i_col<arg_n_cols;i_col++)
	  {
	  Map[i_row][i_col]=NULL;
      } // end for i_col
    } // end for i_row
  // Make space to hold the loaded row details
  n_land_in_row= AllocateLong1DArray(arg_n_rows);
  // Raft is now created
  // Create non NULL grid cells for land
  for(i_row=0;i_row<arg_n_rows;i_row++)
	{
	i_land_count=0;
	for(i_col=0;i_col<arg_n_cols;i_col++)
	  {
	  if(arg_s_condition[i_row][i_col]>0)
		{
		Map[i_row][i_col]=new TGridCell();
		// count valid cells
		i_land_count++;
		} // end if arg_s_condition[i_row][i_col]>0
	  } // end for i_col
	n_land_in_row[i_row]=i_land_count;
	} // end for i_row

 // initialise all the data time series and load into the raft
 //InitialiseAllSeries();

} // end constr 2

//---------------------------------------------------------------------------
   // destructor    // NEED TO REVISE THIS ********
TRaft::~TRaft(void)
{
long int i_row; // row counter
long int i_col;// col counter

  for(i_row=0;i_row<n_rows;i_row++)
    {
    for(i_col=0;i_col<n_cols;i_col++)
      {
      if(Map[i_row][i_col]!=NULL)
        delete (Map[i_row][i_col]);
      }
    }
  // Free up both levels of the allocated array
  delete [] *Map;
  delete [] Map;

  FreeLong1DArray(n_land_in_row);

}

//---------------------------------------------------------------------------
//*****************************************************************************
//* Load Empty Raft Data    // KM May 2014
//*
//* ALLOCATES SPACE FOR DATA REQUIRED TO BE HELD ON THE RAFT
//* Argument; arg_i_row... the row to load
//*
//*****************************************************************************/
void TRaft::LoadEmptyRaftData(int arg_n_rows,
							  int arg_n_cols,
							  int arg_n_env_variables)
{
  // Declare local variables
  long int i_row, i_col, i_env;
  float* f_env_data_array;
  double* d_env_current_z_array;
  double* d_env_current_Vp_array;
  double* d_env_next_z_array;
  double* d_env_next_Vp_array;

  // Loop through the Raft, allocate the necessary grids where required, fill
  // them with zeros
  for(i_row=0;i_row<arg_n_rows;i_row++)
	{
	for(i_col=0;i_col<arg_n_cols;i_col++)
	  {
	  if(Map[i_row][i_col]!=NULL)
		{
		f_env_data_array=AllocateFloat1DArray(arg_n_env_variables);
		Map[i_row][i_col]->SetEnvData(f_env_data_array);
		d_env_current_z_array=AllocateDouble1DArray(arg_n_env_variables);
		Map[i_row][i_col]->SetCurrentZ(d_env_current_z_array);
		d_env_current_Vp_array=AllocateDouble1DArray(arg_n_env_variables);
		Map[i_row][i_col]->SetCurrentVp(d_env_current_Vp_array);
		d_env_next_z_array=AllocateDouble1DArray(arg_n_env_variables);
		Map[i_row][i_col]->SetNextZ(d_env_next_z_array);
		d_env_next_Vp_array=AllocateDouble1DArray(arg_n_env_variables);
		Map[i_row][i_col]->SetNextVp(d_env_next_Vp_array);
		// Fill in the allocated arrays with zeros
		for(i_env=0;i_env<arg_n_env_variables;i_env++)
		  {
		  f_env_data_array[i_env]=0;
		  d_env_current_z_array[i_env]=0;
		  d_env_current_Vp_array[i_env]=0;
		  d_env_next_z_array[i_env]=0;
		  d_env_next_Vp_array[i_env]=0;
		  }
		} // end if Map[i_row][i_col]!=NULL
	  } // end for i_col
	} // end for i_row

} // end LoadEmptyRaftData
//*****************************************************************************
//*****************************************************************************
//* FREE Raft Data    // KM May 2014
//*
//* FREES DATA ALLOCATED TO THE RAFT
//* Argument;
//*
//*****************************************************************************
void TRaft::FreeRaftData(int arg_n_rows,
						 int arg_n_cols)
{
  // Declare local variables
  long int i_row, i_col;
  float* f_env_data_array;
  double* d_env_current_z_array;
  double* d_env_current_Vp_array;
  double* d_env_next_z_array;
  double* d_env_next_Vp_array;

  // Loop through the Raft, allocate the necessary grids where required, fill
  // them with zeros
  for(i_row=0;i_row<arg_n_rows;i_row++)
	{
	for(i_col=0;i_col<arg_n_cols;i_col++)
	  {
	  if(Map[i_row][i_col]!=NULL)
		{
		// The Env data array
		f_env_data_array= Map[i_row][i_col]->GetEnvData();
		if(f_env_data_array!=NULL)
		  FreeFloat1DArray(f_env_data_array);
		f_env_data_array=NULL;
		// The Env current z array
		d_env_current_z_array= Map[i_row][i_col]->GetCurrentZ();
		if(d_env_current_z_array!=NULL)
		  FreeDouble1DArray(d_env_current_z_array);
		d_env_current_z_array=NULL;
		// The Env current Vp array
		d_env_current_Vp_array= Map[i_row][i_col]->GetCurrentVp();
		if(d_env_current_Vp_array!=NULL)
		  FreeDouble1DArray(d_env_current_Vp_array);
		d_env_current_Vp_array=NULL;
		// The Env next z array
		d_env_next_z_array= Map[i_row][i_col]->GetNextZ();
		if(d_env_next_z_array!=NULL)
		  FreeDouble1DArray(d_env_next_z_array);
		d_env_next_z_array=NULL;
		// The Env next Vp array
		d_env_next_Vp_array= Map[i_row][i_col]->GetNextVp();
		if(d_env_next_Vp_array!=NULL)
		  FreeDouble1DArray(d_env_next_Vp_array);
		d_env_next_Vp_array=NULL;
		} // end if Map[i_row][i_col]!=NULL
	  } // end for i_col
	} // end for i_row

} // end FreeRaftData
//*****************************************************************************

//*****************************************************************************
//* Read In Whole ENVStack
//* Loads the supplied file to its position on the raft
//* The file must be in the semi RLE encoded *.MCG format, with RLE for the NODATA
//* Doesn't allocate--just overwrites in each TGridCell
//* Arguments:
//       string arg_s_file  the full path to load the FLOAT MCG file from
//* Returns: ONE success, ZERO Failure
//*****************************************************************************
int TRaft::ReadInWholeENVStack(string arg_s_file,
							   int arg_n_env_variables)
{
int b_success=FALSE;  // return flag
ifstream InFileStream;
int i;
short* i_s_ptr;
short i_s;
long int i_col;
long int i_row;
int i_layer;            // counter
float* f_ptr;   // reading pointer
float f_buffer; // place to read stuff we don't care about.

  i_s_ptr=&(i_s);
  i_col=MINUS_ONE;  // first increment =0
  i_row=0;
  // Open binary file up for binary reading
  InFileStream.open(arg_s_file.c_str(), ios::in |ios::binary);
  while(!InFileStream.eof())
	{
    // Read run length and type
	InFileStream.read((char*)i_s_ptr,sizeof(short));
    if(i_s<0)
      {
	  if(i_s<=K_EOF) // end file
        break;
      // nodata run
      i_s=0-i_s; // make positive
	  i_col+=i_s;

	  if(i_col==n_cols-1)
        {
		i_row++;
		i_col=MINUS_ONE;  // first increment =0
        }
      }
    else
      {
      // A run of data, read out layer by layer
	  for(i=0;i<i_s;i++)
		{
        i_col++;
		if(Map[i_row][i_col]==NULL)
          {
          // NO DATA CELL HERE>>read out the file data to skip this cell
		  for(i_layer=0;i_layer<arg_n_env_variables;i_layer++)
			{
			f_ptr=&f_buffer;
			InFileStream.read((char*)f_ptr,sizeof(float));
            }
          }
        else
		  {
          // Read in all the layers direct to the cell  NO ALLOC
		  Map[i_row][i_col]->SetENVStack(&InFileStream,
										 arg_n_env_variables);
          }
		}  // end for i
	  if(i_col==n_cols-1)
        {
        i_row++;
        i_col=MINUS_ONE; // first increment =0
        }
      } // end elsedata
    if(i_row==n_rows)
      break;
    } // end while eof
  //Close input file
  if(InFileStream.is_open())
    InFileStream.close();

return b_success;
}// end func Read In whole GDMStack


//*****************************************************************************

