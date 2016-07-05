#include <Rcpp.h>
using namespace Rcpp;
#include "AdaptorFileUtils_V2.h"
#include "CompactorFunctions.h"

/////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
int rcpp_AdaptR(std::string PfileName)
  {
  std::string str_file;
  int i_complete = 0;
  str_file = PfileName; // Get the parameter filename from the argument to main
  str_file = str_file.c_str(); // Ensure the parameter file name is a string  ??
  i_complete = Run_Adaptor(str_file); // Call to the function "Run_Adaptor" which does stuff, using the parameter file (str_file) as input
  return i_complete;
  }  // end rcpp_AdaptR

/////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
int rcpp_MuruCompactor(std::string PfileName)
{
  int i_result = 0;
  int i_param = 1; 
  std::string str_file;
  str_file = PfileName; // Get the parameter filename from the argument to main
  str_file = str_file.c_str(); // Ensure the parameter file name is a string  ??   
  i_result = ConvertASCIIGridStackToCompactFloatBinaryFile(i_param,
                                                           str_file);
  return i_result;
} // end rcpp_MuruCompactor

/////////////////////////////////////////////////////////////////////////////////////////////

