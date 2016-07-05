//*****************************************************************************
// Dynamic Array . cpp
// Suite of functions for the allocation and deletion of dynamic arrays, 1& 2D
// Written Tom Harwood May 2010
//*****************************************************************************

#include "DynamicArray.h"

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
// ***** Allocate Double 1D Array
// Creates a 1D array of the specified dimension
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeDouble1DArray(float**)
// Arguments:
// const int arg_i_n_cols    number of columns
// Returns:
// Pointer  to the double* pointing to the newly allocated 1D array
//*****************************************************************************
double* AllocateDouble1DArray(const int arg_i_n_cols)
{
  double* d_new_ptr = new double[arg_i_n_cols];
  return d_new_ptr;
} // End func AllocateDouble1DArray
//*****************************************************************************
// ***** Free Double 1D Array ********
// Deletes the suuplied 1D array of the supplied TYPE
// Memory: Deletes argument
// Arguments:
// T* arg_T_Array    the array to delete  // any type
// Returns: none
//
//*****************************************************************************
void FreeDouble1DArray(double* arg_d_Array)
{
  delete [] arg_d_Array;
}  // end func delete double 1D array
//*****************************************************************************



//*****************************************************************************
//*****************************************************************************
// ***** Allocate Float 1D Array
// Creates a 1D array of the specified dimension
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeFloat1DArray(float**)
// Arguments:
// const int arg_i_n_cols    number of columns
// Returns:
// Pointer  to the float* pointing to the newly allocated 1D array
//*****************************************************************************
float* AllocateFloat1DArray(const int arg_i_n_cols)
{
  float* f_new_ptr = new float[arg_i_n_cols];
  return f_new_ptr;
} // End func AllocateInt1DArray
//*****************************************************************************
// ***** Free INT 1D Array ********
// Deletes the suuplied 1D array of the supplied TYPE
// Memory: Deletes argument
// Arguments:
// T* arg_T_Array    the array to delete  // any type
// Returns: none
//
//*****************************************************************************
void FreeFloat1DArray(float* arg_f_Array)
{
  delete [] arg_f_Array;
}  // end func delete float 1D array
//*****************************************************************************


//*****************************************************************************
//*****************************************************************************
// ***** Allocate Int 1D Array
// Creates a 1D array of the specified dimension
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeInt1DArray(int**)
// Arguments:
// const int arg_i_n_cols    number of columns
// Returns:
// Pointer  to the TEMPLATE type pointing to the newly allocated 1D array
//*****************************************************************************
int* AllocateInt1DArray(const int arg_i_n_cols)
{
  // Allocate memory for array
  int* i_new_ptr = new int[arg_i_n_cols];
  return i_new_ptr;
} // End func AllocateInt1DArray
//*****************************************************************************
// ***** Free INT 1D Array ********
// Deletes the suuplied 1D array of the supplied TYPE
// Memory: Deletes argument
// Arguments:
// T* arg_T_Array    the array to delete  // any type
// Returns: none
//
//*****************************************************************************
void FreeInt1DArray(int* arg_i_Array)
{
  // Free up array
  delete [] arg_i_Array;
}  // end func delete INT 1D array
//*****************************************************************************
//*****************************************************************************
// ***** Allocate Long 1D Array
// Creates a 1D array of the specified dimension
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeInt1DArray(int**)
// Arguments:
// const long int arg_i_n_cols    number of columns
// Returns:
// Pointer  to the TEMPLATE type pointing to the newly allocated 1D array
//*****************************************************************************
long int* AllocateLong1DArray(const long int arg_i_n_cols)
{
  // Allocate memory for array
  long int* i_new_ptr = new long int[arg_i_n_cols];
  return i_new_ptr;
} // End func AllocateLong1DArray
//*****************************************************************************
// ***** Free Long 1D Array ********
// Deletes the suuplied 1D array of the supplied TYPE
// Memory: Deletes argument
// Arguments:
// T* arg_T_Array    the array to delete  // any type
// Returns: none
//
//*****************************************************************************
void FreeLong1DArray(long int* arg_i_Array)
{
  // Free up array
  delete [] arg_i_Array;
}  // end func delete LONG INT 1D array
//*****************************************************************************


//*****************************************************************************
//*****************************************************************************
// ***** Allocate Short 1D Array
// Creates a 1D array of the specified dimension
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeShort1DArray(short**)
// Arguments:
// const short arg_i_n_cols    number of columns
// Returns:
// Pointer  to the TEMPLATE type pointing to the newly allocated 1D array
//*****************************************************************************
short* AllocateShort1DArray(const int arg_i_n_cols)
{
  // Allocate memory for array
  short* s_new_ptr = new short[arg_i_n_cols];
  return s_new_ptr;
} // End func AllocateShort1DArray
//*****************************************************************************
// ***** Free SHORT 1D Array ********
// Deletes the suuplied 1D array of the supplied TYPE
// Memory: Deletes argument
// Arguments:
// T* arg_T_Array    the array to delete  // any type
// Returns: none
//
//*****************************************************************************
void FreeShort1DArray(short* arg_s_Array)
{
  // Free up array
  delete [] arg_s_Array;
}  // end func delete Short 1D array
//*****************************************************************************

//*****************************************************************************
//*****************************************************************************
// ***** Allocate Double 2D Array ********
// Creates a 2D array of the specified dimensions
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeDouble2DArray(double**)
// Arguments:
// const int arg_i_n_rows    number of rows
// const int arg_i_n_cols    number of columns
// Returns:
// Pointer to Pointer to Float pointing to the newly allocated 2D array
// ..which can be referenced as result[][].
// usage example:
//      double** d_fish;
//      d_fish= AllocateDouble2DArray(10,10);
//      d_fish[2][5]=21.563;
//      FreeDouble2DArray(d_fish);
//*****************************************************************************
double **AllocateDouble2DArray(const int arg_i_n_rows,
                               const int arg_i_n_cols)
{
int i_row;              // row counter
double** d_new_pptr;     // Pointer to pointer for the whole 2D buffer
double* d_row_ptr;       // Pointer to start of each row in the array

  // Allocate memory for array of elements of column ..
  // i.e. list of pointers to first element of each row
  d_new_pptr = new double*[arg_i_n_rows];
  // Allocate memory for array of elements of all the rows
  // Do this as one big block to speed things up
  d_row_ptr = new double [arg_i_n_rows* arg_i_n_cols];

  // Now point the pointers in the right place
  for(i_row=0; i_row<arg_i_n_rows; i_row++)
    {
    *(d_new_pptr+i_row) = d_row_ptr;
    d_row_ptr += arg_i_n_cols;
    }
  return d_new_pptr;
}  // end func Allocate Double 2D Array

//*****************************************************************************
// ***** Free Double 2D Array ********
// Deletes the suuplied 2D array
// Memory: Deletes argument
// Arguments:
// double** arg_d_Array    the array to delete
// Returns: none
//
//*****************************************************************************
void FreeDouble2DArray(double** arg_d_Array)
{
  // Free up both levels of the allocated array
  delete [] *arg_d_Array;
  delete [] arg_d_Array;
} // end func delete Double 2D array

//*****************************************************************************
//*****************************************************************************
// ***** Allocate Float 2D Array ********
// Creates a 2D array of the specified dimensions
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeFloat2DArray(float**)
// Arguments:
// const int arg_i_n_rows    number of rows
// const int arg_i_n_cols    number of columns
// Returns:
// Pointer to Pointer to Float pointing to the newly allocated 2D array
// ..which can be referenced as result[][].
// usage example:
//      float** f_fish;
//      f_fish= AllocateFloat2DArray(10,10);
//      f_fish[2][5]=21.563;
//      FreeFloat2DArray(f_fish);
//*****************************************************************************
float **AllocateFloat2DArray( const int arg_i_n_rows,
                              const int arg_i_n_cols)
{
int i_row;              // row counter
float** f_new_pptr;     // Pointer to pointer for the whole 2D buffer
float* f_row_ptr;       // Pointer to start of each row in the array

  // Allocate memory for array of elements of column ..
  // i.e. list of pointers to first element of each row
  f_new_pptr = new float*[arg_i_n_rows];
  // Allocate memory for array of elements of all the rows
  // Do this as one big block to speed things up
  f_row_ptr = new float [arg_i_n_rows* arg_i_n_cols];

  // Now point the pointers in the right place
  for(i_row=0; i_row<arg_i_n_rows; i_row++)
    {
    *(f_new_pptr+i_row) = f_row_ptr;
    f_row_ptr += arg_i_n_cols;
    }
  return f_new_pptr;
}  // end func Allocate Float 2D Array

//*****************************************************************************
// ***** Free Float 2D Array ********
// Deletes the suuplied 2D array
// Memory: Deletes argument
// Arguments:
// float** arg_f_Array    the array to delete
// Returns: none
//
//*****************************************************************************
void FreeFloat2DArray(float** arg_f_Array)
{
  // Free up both levels of the allocated array
  delete [] *arg_f_Array;
  delete [] arg_f_Array;
} // end func delete Float 2D array



//*****************************************************************************
//*****************************************************************************
// ***** Allocate Int 2D Array ********
// Creates a 2D array integers of the specified dimensions
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeInt2DArray(int**)
// Arguments:
// const int arg_i_n_rows    number of rows
// const int arg_i_n_cols    number of columns
// Returns:
// Pointer to Pointer to Int pointing to the newly allocated 2D array
// ..which can be referenced as result[][].
// usage example:
//      int** i_fish;
//      i_fish= AllocateInt2DArray(60,300);
//      i_fish[27][35]=26.563;
//      FreeInt2DArray(i_fish);
//*****************************************************************************
int** AllocateInt2DArray( const int arg_i_n_rows,
                            const int arg_i_n_cols)
{
int i_row;              // row counter
int** i_new_pptr;     // Pointer to pointer for the whole 2D buffer
int* i_row_ptr;       // Pointer to start of each row in the array

  // Allocate memory for array of elements of column ..
  // i.e. list of pointers to first element of each row
  i_new_pptr = new int*[arg_i_n_rows];
  // Allocate memory for array of elements of all the rows
  // Do this as one big block to speed things up
  i_row_ptr = new int[arg_i_n_rows*arg_i_n_cols];

  // Now point the pointers in the right place
  for(i_row=0; i_row<arg_i_n_rows; i_row++)
    {
    *(i_new_pptr+i_row) = i_row_ptr;
    i_row_ptr += arg_i_n_cols;
    }
  return i_new_pptr;
}  // end func Allocate Int 2D Array
//*****************************************************************************
// ***** Free Int 2D Array ********
// Deletes the suuplied 2D array
// Memory: Deletes argument
// Arguments:
// int** arg_i_Array    the array to delete
// Returns: none
//
//*****************************************************************************
void FreeInt2DArray(int** arg_i_Array)
{
  // Free up both levels of the allocated array
  delete [] *arg_i_Array;
  delete [] arg_i_Array;
} // end func delete Int 2D array






//*****************************************************************************
// ***** Allocate Short 2D Array ********
// Creates a 2D array short integers of the specified dimensions
// Memory: ###### Allocates an array which must be DELETED ELSEWHERE! #######
// use FreeShort2DArray(short**)
// Arguments:
// const short arg_i_n_rows    number of rows
// const short arg_i_n_cols    number of columns
// Returns:
// Pointer to Pointer to short pointing to the newly allocated 2D array
// ..which can be referenced as result[][].
// usage example:
//      short** s_fish;
//      s_fish= AllocateShort2DArray(60,300);
//      s_fish[27][35]=26.563;
//      FreeShort2DArray(s_fish);
//*****************************************************************************
short** AllocateShort2DArray(const int arg_i_n_rows,
                             const int arg_i_n_cols)
{
int i_row;              // row counter
short** s_new_pptr;     // Pointer to pointer for the whole 2D buffer
short* s_row_ptr;       // Pointer to start of each row in the array

  // Allocate memory for array of elements of column ..
  // i.e. list of pointers to first element of each row
  s_new_pptr = new short*[arg_i_n_rows];
  // Allocate memory for array of elements of all the rows
  // Do this as one big block to speed things up
  s_row_ptr = new short[arg_i_n_rows*arg_i_n_cols];

  // Now point the pointers in the right place
  for(i_row=0; i_row<arg_i_n_rows; i_row++)
    {
    *(s_new_pptr+i_row) = s_row_ptr;
    s_row_ptr += arg_i_n_cols;
    }
  return s_new_pptr;
}  // end func Allocate Short 2D Array
//*****************************************************************************
// ***** Free short 2D Array ********
// Deletes the suuplied 2D array
// Memory: Deletes argument
// Arguments:
// short** arg_s_Array    the array to delete
// Returns: none
//
//*****************************************************************************
void FreeShort2DArray(short** arg_s_Array)
{
  // Free up both levels of the allocated array
  delete [] *arg_s_Array;
  delete [] arg_s_Array;
} // end func delete short 2D array



//---------------------------------------------------------------------------

 