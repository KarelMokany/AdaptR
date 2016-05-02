#' Run AdaptR
#'
#' @param Pfile A file path to the parameter file
#' @return A set of output files in the specified location
#' @examples
#' CompactGrids("C:/MyDocuments/AdaptR_ParameterFile.txt")
AdaptR<-function(Pfile = Parameter.File.Path)
{
  # load the dll
  dyn.load("src/main.dll")
  # call AdaptR dll
  .C("AdaptR",  argv = as.character(c(Pfile)),
     arg_i_catch = as.integer(c(0,0)))
} # end AdaptR function



