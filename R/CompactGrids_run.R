
#' Create compacted environmental grids for AdaptR
#'
#' @param Pfile A file path to the parameter file
#' @return A set of cpmacted environment files in the specified location
#' @examples
#' CompactGrids("C:/MyDocuments/CompactorParameterFile.txt")
CompactGrids<-function(Pfile = Parameter.File.Path)
{
  # load the dll
  dyn.load("src/main.dll")
  # call MuruCompactor from dll
  .C("MuruCompactor",  argv = as.character(c(Pfile)),
     arg_i_catch = as.integer(c(0,0)))
  unload("src/main.dll")
} # end CompactGrids function