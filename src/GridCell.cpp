//*****************************************************************************
// TGridCell . CPP:     (MURU)
// written by Tom Harwood Jan-Aug 2011
// base class for data which is stored for land.. (NULL for sea)
// Handles multiple time steps via juggled pairs of stacks, some of which are
// dynamically allocated.
// Use of increment to interpolate between known data steps
// Can take allocated species data
//
//*****************************************************************************
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include "defines.h"
//#include "structures.h"
#include "DynamicArray.h"
#include "GridCell.h"

using namespace std;
//****************************************************************************
// Constructor 1 initialises things
//****************************************************************************
TGridCell::TGridCell(void)
{

  // Nullify pointers
  f_env_data=NULL;  		// pointer to a dynamically allocated array to hold the values of each environmental variable for this cell
  d_env_current_z=NULL; // pointer to a dynamically allocated array to hold the values of the mean tolerance attribute (z) for each environment variable
  d_env_current_Vp=NULL; // pointer to a dynamically allocated array to hold the values of the variance in tolerance attribute (Vp) for each environment variable
  d_env_next_z=NULL;
  d_env_next_Vp=NULL;

  // Zero the other holders
  s_time_last_present=0; // holds the last timestep this species was present in this grid cell
  s_present=0;
  s_newly_colonised=0;

} // end constr 1
/******************************************************************************
 * TGridCell  DESTRUCTOR
 * deletes stuff
 *********************************************************************/
TGridCell::~TGridCell( void)
{
  if(f_env_data!=NULL)
	delete [] f_env_data;
  if(d_env_current_z!=NULL)
	delete [] d_env_current_z;
  if(d_env_current_Vp!=NULL)
	delete [] d_env_current_Vp;
  if(d_env_current_z!=NULL)
	delete [] d_env_next_z;
  if(d_env_current_Vp!=NULL)
	delete [] d_env_next_Vp;

}              // end destructor

/******************************************************************************
 * TGridCell  Set ENVStack
 * loads data for this grid cell from the supplied file stream
 * Arguments
 * const int arg_i_GDMStack,   //K_FROM_STACK 0 and K_TO_STACK 1
 * ifstream arg_InFileStream,  // open binary file stream to read FLOATS from
 * const int arg_n_layers      // number of layers
 * We require the user of the function to pass in a file stream which is a
 * FLOAT grid by reading the header file first....
 *****************************************************************************/
void TGridCell::SetENVStack(ifstream* arg_InFileStream_ptr,
                            const int arg_n_layers)
{
int i_layer;            // counter
float* f_ptr;   // reading pointer
//float f;
  for(i_layer=0;i_layer<arg_n_layers;i_layer++)
	{
	f_ptr=&f_env_data[i_layer];
	(*arg_InFileStream_ptr).read((char*)f_ptr,sizeof(float));
	}
} // end func Set GDMStack
//******************************************************************************
