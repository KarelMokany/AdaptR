//*****************************************************************************
// Dynamic Array . h
// Suite of functions for the allocation and deletion of dynamic arrays, 1& 2D
// Written Tom Harwood May 2010
// Shorts added by K Mokany May 2010
//*****************************************************************************

#ifndef DynamicArrayH
#define DynamicArrayH
//*****************************************************************************
//*****************************************************************************
//***1 D Arrays*******************
// Double
double* AllocateDouble1DArray(const int arg_i_n_cols);
void FreeDouble1DArray(double* arg_d_Array);
// Float
float* AllocateFloat1DArray(const int arg_i_n_cols);
void FreeFloat1DArray(float* arg_f_Array);
// Integer
int* AllocateInt1DArray(const int arg_i_n_cols);
void FreeInt1DArray(int* arg_i_Array);
// Long int
long int* AllocateLong1DArray(const long int arg_i_n_cols);
void FreeLong1DArray(long int* arg_i_Array);
// Short
short* AllocateShort1DArray(const int arg_i_n_cols);
void FreeShort1DArray(short* arg_s_Array);

//***2 D Arrays*******************
// Double
double** AllocateDouble2DArray(const int arg_i_n_rows,
                               const int arg_i_n_cols);
void FreeDouble2DArray(double** arg_d_Array);

// Float
float **AllocateFloat2DArray( const int arg_i_n_rows,
                              const int arg_i_n_cols);
void FreeFloat2DArray(float** arg_f_Array);
// Integer
int** AllocateInt2DArray( const int arg_i_n_rows,
                          const int arg_i_n_cols);
void FreeInt2DArray(int** arg_i_Array);
// Short
short** AllocateShort2DArray(const int arg_i_n_rows,
                             const int arg_i_n_cols);
void FreeShort2DArray(short** arg_s_Array);


//---------------------------------------------------------------------------
#endif
 