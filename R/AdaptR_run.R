################################################################################
##
## Run AdaptR - Preliminary code
##
## Karel Mokany - 28 April 2016
##
################################################################################

# Set the working directory
#setwd("C:/Users/mok010/Karels Files/Mokany Files/My R applications/AdaptR/source")



## Set the parameter file name
#Parameter.File.Path = "C:/Users/mok010/Karels Files/Mokany Files/My R applications/AdaptR/Test_Inputs/Parameters/Adaptor_parameters_ironensis_TEST.txt"

## Run AdaptR
AdaptR<-function(Pfile = Parameter.File.Path)
{
  # load the dll
  dyn.load("src/main.dll")
  # call AdaptR dll
  .C("AdaptR",  argv = as.character(c(Pfile)),
     arg_i_catch = as.integer(c(0,0)))
} # end AdaptR function

## Compact Environment Files
CompactGrids<-function(Pfile = Parameter.File.Path)
{
  # load the dll
  dyn.load("src/main.dll")
  # call MuruCompactor from dll
  .C("MuruCompactor",  argv = as.character(c(Pfile)),
     arg_i_catch = as.integer(c(0,0)))
} # end CompactGrids function

