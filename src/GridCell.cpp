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
 * TGridCell  Set Freeze
 * Freezes and melts the FROM stack by storing a local copy of outdated data
 * Set arg_f_freeze to K_MELT to unfreeze or K_FREEZE to make a copy of the
 * current stack
 * Arguments: const float arg_f_freeze,      // K_MELT to unfreeze or K_FREEZE
                          const int arg_i_from_stack, // needed if FREEZING
                          const int arg_n_layers      // needed if FREEZING
 *********************************************************************
void TGridCell::SetFreeze(const float arg_f_freeze,
                          const int arg_i_from_stack,
                          const int arg_n_layers)
{
int i_layer;
  if(arg_f_freeze==K_MELT)
    {
    f_gdm_data_ptr[K_FREEZE][0]=K_MELT;
    }
  else
    {
    // Freeze to the current FROM layer
    for(i_layer=0;i_layer<arg_n_layers;i_layer++)
      {
      f_gdm_data_ptr[K_FREEZE][i_layer]= f_gdm_data_ptr[arg_i_from_stack][i_layer];
      } // end for i_layer
    }
}  // end func Set Freeze
/******************************************************************************
 * TGridCell  Set GDMStack
 * loads data for this grid cell from the supplied file stream
 * Arguments
 * const int arg_i_GDMStack,   //K_FROM_STACK 0 and K_TO_STACK 1
 * ifstream arg_InFileStream,  // open binary file stream to read FLOATS from
 * const int arg_n_layers      // number of layers
 * We require the user of the function to pass in a file stream which is a
 * FLOAT grid by reading the header file first....
 *********************************************************************
void TGridCell::SetGDMStack(const int arg_i_GDMStack,
                            ifstream* arg_InFileStream_ptr,
                            const int arg_n_layers)
{
int i_layer;            // counter
float* f_ptr;   // reading pointer
//float f;
  for(i_layer=0;i_layer<arg_n_layers;i_layer++)
    {
    f_ptr=&f_gdm_data_ptr[arg_i_GDMStack][i_layer];
    (*arg_InFileStream_ptr).read((char*)f_ptr,sizeof(float));
    }
} // end func Set GDMStack
/******************************************************************************
 * TGridCell  Allocate GDMStack
 * Allocates a stack of n_layers at f_gdm_data_ptr[arg_i_GDMStack]
 *********************************************************************
void TGridCell::AllocateGDMStack(const int arg_i_GDMStack,
                                 const int arg_n_layers)  // allocates space for the stack at this cell
{
  // Allocate a 1 D array of floats at f_gdm_data_ptr[arg_i_GDMStack]
  f_gdm_data_ptr[arg_i_GDMStack]=AllocateFloat1DArray(arg_n_layers);
} // end func Allocate GDMStack
/******************************************************************************
 * TGridCell  Free GDMStack
 * deletes the indexed memory stack
 *********************************************************************
void TGridCell::FreeGDMStack(const int arg_i_GDMStack)  // frees space for the stack at this cell
{
  if(f_gdm_data_ptr[arg_i_GDMStack]!=NULL)
    FreeFloat1DArray(f_gdm_data_ptr[arg_i_GDMStack]);
  f_gdm_data_ptr[arg_i_GDMStack]=NULL;
}// end func free GDMStack

/******************************************************************************
 * TGridCell  Load GDMStack
 * Frees up any allocation, allocates and then loads in the stack of numbers
 * Arguments
 * const int arg_i_GDMStack,   //K_FROM_STACK 0 and K_TO_STACK 1
 * ifstream arg_InFileStream,  // open binary file stream to read from
 * const int arg_n_layers      // number of layers allocation size
 *********************************************************************
void TGridCell::LoadGDMStack( const int arg_i_GDMStack,   // Allocates and sets
                              ifstream* arg_InFileStream_ptr,
                              const int arg_n_layers)
{
  FreeGDMStack(arg_i_GDMStack);     // frees up GDMStack if allocated
  AllocateGDMStack(arg_i_GDMStack,
                   arg_n_layers); // Allocates memory
  SetGDMStack(arg_i_GDMStack,       //Assume all is well and just set
              arg_InFileStream_ptr,
              arg_n_layers);
}// end func load GDMStack




/******************************************************************************
 * TGridCell  Allocate Incrememt
 * Allocates a stack of n_layers at f_increment
 *********************************************************************
void TGridCell::AllocateGDMIncrement(const int arg_n_layers)  // allocates space for the stack at this cell
{
  f_gdm_increment=AllocateFloat1DArray(arg_n_layers);
} //end func Allocate Increment
/******************************************************************************
 * TGridCell  Free Incrememt
 * Frees a stack of n_layers at f_increment
 *********************************************************************
void TGridCell::FreeGDMIncrement(void)
{
  if(f_gdm_increment!=NULL)
    FreeFloat1DArray(f_gdm_increment);
  f_gdm_increment=NULL;
} //end func Free Increment

/******************************************************************************
 * TGridCell  Calculate GDM Increment
 * calculates the increment for each layer in the stacks
 * stores results in f_increment.
 * While we are doing this, we RESET the value of the TO stack to be FROM_STACK + INCREMENT
 * NOTE f_increment is ALLOCATED so may not be there
 *  Arguments:
 * const int arg_i_step,  // number of steps between layers.. e.g. 2020-2030 arg_i_step=10
 * const int arg_n_layers         // number of GDM layers
 * const int arg_i_from_stack // so increment goes in the right direction!
 *  Returns:
 *   1   OK  everything OK
 *   0   BAD f_increment not allocated
 *********************************************************************
int TGridCell::CalcGDMIncrement(const int arg_i_step,   // Calculates arg_i_step increment between each stack
                                const int arg_n_layers,
                                const int arg_i_from_stack) // number of layers
{
int i_layer;     // loop counter
int i_to_stack;
  // Test for allocation
  if(f_gdm_increment==NULL)
    return BAD;

    if(arg_i_from_stack==0)
      i_to_stack=1;
    else
      i_to_stack=0;
  // Loop through all layers calculating increment and saving in f_increment
  // NOTE the from and to stacks can leap around, so we must take a value as an argument
  for(i_layer=0;i_layer<arg_n_layers;i_layer++)
    {
    if(arg_i_step>0)
      f_gdm_increment[i_layer]=(f_gdm_data_ptr[i_to_stack][i_layer]-f_gdm_data_ptr[arg_i_from_stack][i_layer])/(float)arg_i_step;
    else
      f_gdm_increment[i_layer]=0;
    f_gdm_data_ptr[i_to_stack][i_layer]=f_gdm_data_ptr[arg_i_from_stack][i_layer]+f_gdm_increment[i_layer];
    }
  return OK;

} // end func CalcIncrement
/******************************************************************************
 * TGridCell  Increment GDM Increment
 * adds on the increment for each layer to the TO stack
 * stores results in f_increment.
 * WE HAVE JUST SWAPPED USING FLIPSTACKS (ideally) so
 * the from stack has last years TO stack in it
 * Add increment to the FROM stack and store in the TO stack
 * NOTE f_increment is ALLOCATED so may not be there
 *  Arguments:
 * const int arg_i_to_stack,  // stack to add to
 *********************************************************************
void TGridCell::IncrementGDMIncrement(const int arg_i_to_stack,
                                      const int arg_n_layers) // number of layers
{
int i_layer;     // loop counter
int i_from_stack;

  if(arg_i_to_stack>0)
    i_from_stack=0;
  else
    i_from_stack=1;

  // Loop through all layers adding increment and saving in f_increment
  for(i_layer=0;i_layer<arg_n_layers;i_layer++)
    f_gdm_data_ptr[arg_i_to_stack][i_layer]=f_gdm_data_ptr[i_from_stack][i_layer]+f_gdm_increment[i_layer];
} // end func Increment GDM increment

/******************************************************************************
 * TGridCell  Calculate RICH Increment
 * calculates the increment for each layer in the stacks
 * stores results in f_increment.
 * While we are doing this, we RESET the value of the TO stack to be FROM_STACK + INCREMENT
 * NOTE f_increment is ALLOCATED so may not be there
 *  Arguments:
 * const int arg_i_step,  // number of steps between layers.. e.g. 2020-2030 arg_i_step=10
 * const int arg_i_from_stack // so increment goes in the right direction!
 *  Returns:
 *   1   OK  everything OK
 *   0   BAD f_increment not allocated
 *********************************************************************
int TGridCell::CalcRichIncrement(const int arg_i_step,   // Calculates arg_i_step increment between each stack
                                const int arg_i_from_stack) // current from stack
{
int i_to_stack;

    if(arg_i_from_stack==0)
      i_to_stack=1;
    else
      i_to_stack=0;
  // Calculate increment and saving in f_increment
  // NOTE the from and to stacks can leap around, so we must take a value as an argument
  if(arg_i_step>0)
    {
    f_rich_increment=(f_richness[i_to_stack]-f_richness[arg_i_from_stack])/(float)arg_i_step;
    if(fabs(f_rich_increment)<NEARLY_ZERO)
      f_rich_increment=0;
    }
  else
    f_rich_increment=0;
  f_richness[i_to_stack]=f_richness[arg_i_from_stack]+f_rich_increment;

  return OK;

} // end func CalcRichIncrement
/******************************************************************************
 * TGridCell  Increment Rich Increment
 * adds on the increment to the TO stack
 *  Arguments:
 * const int arg_i_to_stack,  // stack to add to
 *********************************************************************
void TGridCell::IncrementRichIncrement(const int arg_i_to_stack)
{
int i_from_stack;

  if(arg_i_to_stack>0)
    i_from_stack=0;
  else
    i_from_stack=1;

    f_richness[arg_i_to_stack]=f_richness[i_from_stack]+f_rich_increment;
} // end func Increment Rich increment


/******************************************************************************
 * TGridCell  Calculate Cond Increment
 * calculates the increment stores results in f_cond_increment.
 * While we are doing this, we RESET the value of the TO stack to be FROM_STACK + INCREMENT
 *  Arguments:
 * const int arg_i_step,  // number of steps between layers.. e.g. 2020-2030 arg_i_step=10
 * const int arg_i_from_stack // so increment goes in the right direction!
 *  Returns:
 *   1   OK  everything OK
 *   0   BAD f_increment not allocated
 *********************************************************************
int TGridCell::CalcCondIncrement(const int arg_i_step,   // Calculates arg_i_step increment between each stack
                                const int arg_i_from_stack) // current from stack
{
int i_to_stack;

    if(arg_i_from_stack==0)
      i_to_stack=1;
    else
      i_to_stack=0;
  // Calculate increment and saving in f_increment
  // NOTE the from and to stacks can leap around, so we must take a value as an argument
  if(arg_i_step>0)
    {
    f_cond_increment=(f_condition[i_to_stack]-f_condition[arg_i_from_stack])/(float)arg_i_step;
    if(fabs(f_cond_increment)<NEARLY_ZERO)
      f_cond_increment=0;
    }
  else
    f_cond_increment=0;

  f_condition[i_to_stack]=f_condition[arg_i_from_stack]+f_cond_increment;

  return OK;

} // end func CalcCondIncrement
/******************************************************************************
 * TGridCell  Increment Cond Increment
 * adds on the increment to the TO stack
 *  Arguments:
 * const int arg_i_to_stack,  // stack to add to
 *********************************************************************
void TGridCell::IncrementCondIncrement(const int arg_i_to_stack)
{
int i_from_stack;

  if(arg_i_to_stack>0)
    i_from_stack=0;
  else
    i_from_stack=1;

  f_condition[arg_i_to_stack]=f_condition[i_from_stack]+f_cond_increment;
} // end func Increment Cond increment

/******************************************************************************
 * TGridCell  Calculate Similarity  (NO ECOL DISTANCE.. a bit quicker)
 * calculates the similarity between the SELECTED data at this cell with the
 * supplied stack of OTHER data  ( this can actually be any stack)
 * Similarity is calculated as the exponential of MINUS the sum absolute
 * difference between all layers   (e^-sumdiff)
 *  Arguments:
 * float* arg_f_OtherData,  // pointer to 1D array of GDm layer data
 * const int arg_i_stack    // internal stack to compare with
 * const int arg_n_layers         // number of GDM layers
 *********************************************************************
float TGridCell::CalcSimilarity(float* arg_f_OtherData,
                                const int arg_i_stack, //which stack we are
                                const int arg_n_layers,       // how many
                                const float arg_f_intercept)
{
int i_layer;     // loop counter
float f_ecol=0;  // zero the sum here
float f_similarity;

  // Loop through all layers summing ecol distance
  for(i_layer=0;i_layer<arg_n_layers;i_layer++)
    f_ecol+=fabs(arg_f_OtherData[i_layer]-f_gdm_data_ptr[arg_i_stack][i_layer]);

  // we have sum difference.. exp scale it to get similarity and return
  f_similarity=exp(-1*(arg_f_intercept+f_ecol));
  
  return (f_similarity);
}// end func Calc Similarity

/******************************************************************************
 * TGridCell  Calculate Similarity (WITH ECOL DIST.. a bit slower)
 * calculates the similarity between the SELECTED data at this cell with the
 * supplied stack of OTHER data  ( this can actually be any stack)
 * Similarity is calculated as the exponential of MINUS the sum absolute
 * difference between all layers   (e^-sumdiff)
 *  Arguments:
 * float* arg_f_OtherData,  // pointer to 1D array of GDm layer data
 * const int arg_i_stack    // internal stack to compare with
 * const int arg_n_layers         // number of GDM layers
 * float* arg_f_ecol_dist   // pointer so we can get ecol dist back when required:
 *  // set to NULL and nothing will happen but better to use the sister function with no argument
 *********************************************************************
float TGridCell::CalcSimilarity(float* arg_f_OtherData,
                                     const int arg_i_stack, //which stack we are
                                     const int arg_n_layers,
                                     float* arg_f_ecol_dist,   // pointer to float which we can set to the ecol distance
                                      const float arg_f_intercept)
{
int i_layer;     // loop counter
float f_ecol=0;  // zero the sum here
float f_similarity;

  // Loop through all layers summing ecol distance
  for(i_layer=0;i_layer<arg_n_layers;i_layer++)
    f_ecol+=fabs(arg_f_OtherData[i_layer]-f_gdm_data_ptr[arg_i_stack][i_layer]);

  // where required, return the ecological distance
  if(arg_f_ecol_dist!=NULL)
    *arg_f_ecol_dist=f_ecol;

  // we have sum difference.. exp scale it to get similarity  and return
  f_similarity=exp(-1*(arg_f_intercept+f_ecol));

  return (f_similarity);
}// end func Calc Similarity

/******************************************************************************
 * TGridCell  Calculate Dissimilarity  (NO ECOL DISTANCE.. a bit quicker)
 * calculates the dissimilarity between the SELECTED data at this cell with the
 * supplied stack of OTHER data  ( this can actually be any stack)
 * Dissimilarity is calculated as the 1-exponential of MINUS the sum absolute
 * difference between all layers   (e^-sumdiff)
 *  Arguments:
 * float* arg_f_OtherData,  // pointer to 1D array of GDm layer data
 * const int arg_i_stack    // internal stack to compare with
 * const int arg_n_layers         // number of GDM layers
 *********************************************************************
float TGridCell::CalcDissimilarity(float* arg_f_OtherData,
                                     const int arg_i_stack, //which stack we are
                                     const int arg_n_layers,       // how many
                                  const float arg_f_intercept)
{
int i_layer;     // loop counter
float f_ecol=0;  // zero the sum here
float f_similarity;

  // Loop through all layers summing ecol distance
  for(i_layer=0;i_layer<arg_n_layers;i_layer++)
    f_ecol+=fabs(arg_f_OtherData[i_layer]-f_gdm_data_ptr[arg_i_stack][i_layer]);

  // we have sum difference.. exp scale it to get similarity and return
  f_similarity=exp(-1*(arg_f_intercept+f_ecol));

  return (1-f_similarity);
}// end func Calc Dissimilarity

/******************************************************************************
 * TGridCell  Calculate Dissimilarity (WITH ECOL DIST.. a bit slower)
 * calculates the dissimilarity between the SELECTED data at this cell with the
 * supplied stack of OTHER data  ( this can actually be any stack)
 * Dissimilarity is calculated as the 1- exponential of MINUS the sum absolute
 * difference between all layers   (e^-sumdiff)
 *  Arguments:
 * float* arg_f_OtherData,  // pointer to 1D array of GDm layer data
 * const int arg_i_stack    // internal stack to compare with
 * const int arg_n_layers         // number of GDM layers
 * float* arg_f_ecol_dist   // pointer so we can get ecol dist back when required:
 *  // set to NULL and nothing will happen but better to use the sister function with no argument
 *********************************************************************
float TGridCell::CalcDissimilarity(float* arg_f_OtherData,
                                   const int arg_i_stack, //which stack we are
                                     const int arg_n_layers,
                                     float* arg_f_ecol_dist,   // pointer to float which we can set to the ecol distance
                                   const float arg_f_intercept)
{
int i_layer;     // loop counter
float f_ecol=0;  // zero the sum here
float f_similarity;

  // Loop through all layers summing ecol distance
  for(i_layer=0;i_layer<arg_n_layers;i_layer++)
    f_ecol+=fabs(arg_f_OtherData[i_layer]-f_gdm_data_ptr[arg_i_stack][i_layer]);

  // where required, return the ecological distance
  if(arg_f_ecol_dist!=NULL)
    *arg_f_ecol_dist=f_ecol;

  // we have sum difference.. exp scale it to get similarity  and return
  f_similarity=exp(-1*(arg_f_intercept+f_ecol));

  return (1-f_similarity);
}// end func Calc Dissimilarity
//---------------------------------------------------------------------------
*/
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
