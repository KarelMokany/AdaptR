#' Run the AdaptR model
#'
#' @param run.name A suitable name for the simulation run
#' @param parameter.file.name A file path to the parameter file, which will be created & called
#' @param ncols The number of columns on the spatial grid
#' @param nrows The number of rows on the spatial grid
#' @param output.folder.path A valid path to a folder where the results will be written
#' @param verbose.outputs If TRUE, predictions will be written out at each time step
#' @param n.time.points The number of time points in the simulation, which will equal the number of compacted environment grids
#' @param n.env.vars The number of environment variables being considered
#' @param env.vars.names A list of the names of the environment variables being considered
#' @param env.grids.folder.path A valid path to the folder containing the compacted environment grids
#' @param env.grids.name.file A file path to the file holding the names of the compacted environment grids, in order of their implementation in the simulation
#' @param species.initial.grid A file path to the ascii file for the initial occurrence of the focal species, where 1 = present, 0 = absent, -9999 = not valid potential habitat.
#' @param minimum.survival.percentage The percent of current population size below which local extinction will occur following a selection event. The complement of this value is used as the threshold below which evolutionary adaptation will occur.
#' @param resident.population.weighting The percent weighting of a resident population when considering admixture through dispersal from neighbouring populations (weighted by dispersal probability)
#' @param dispersal.neighbourhood.file A file path to the dispersal neighbourhood for AdaptR (.dna) file specifying the dispersal probability in the relative neighbourhood around each cell. See "Create_Dispersal"
#' @param species.location.file A file path to a text file specifying the x and y coordinates for populations (grid cells) where outputs will be written for each generation. Format is space delimited text (.txt) of two columns, with row 1 = n_locations n, row 2 = Longitude Latitude, then n rows specifying the x & y locations to track
#' @param env.lower.thresholds A list of length n.env.vars, specifying the lower tolerance threshold value for each environmental variable
#' @param env.upper.thresholds A list of length n.env.vars, specifying the upper tolerance threshold value for each environmental variable
#' @param env.low.adaptation A logical list (TRUE, FALSE) of length n.env.vars, specifying for each environmental variable whether adaptation on the lower tolerance threshold is being considered 
#' @param env.high.adaptation A logical list (TRUE, FALSE) of length n.env.vars, specifying for each environmental variable whether adaptation on the upper tolerance threshold is being considered
#' @param adapt.limit A list of length n.env.vars specifying the value beyond which the threshold tolerance cannot adapt further 
#' @param heritability A list of length n.env.vars specifying the heritability of threshold tolerance attribute for each environmental variable
#' @param fitness.cost A list of length n.env.vars specifying the 'fitness cost' of adaptation away from the original threshold tolerance values. This is the slope of the linear relation between the intensity of stabilising selection and the squared distance from the original threshold tolerance trait value.
#' @param adapt.threshold.grids A logical list (TRUE, FALSE) of length n.env.vars, specifying for each environmental variable whether spatial grids will be provided specifying the initial threshold tolerance in each grid cell
#' @param adapt.threshold.grid.names Where relevant, a list of length n.env.vars specifying the file paths to the ascii files of the initial threshold tolerance value across the grid
#' @param phenotypic.sd.grid A logical list (TRUE, FALSE) of length n.env.vars, specifying for each environmental variable whether spatial grids will be provided specifying the initial phenotypic standard deviation for threhold tolerances in each grid cell
#' @param phenotypic.sd.grid.names Where relevant, a list of length n.env.vars specifying the file paths to the ascii files of the initial threshold tolerance phenotypic standard deviation value across the grid
#' @param phenotypic.sd.value Where relevant, a list of length n.env.vars specifying the initial phenotypic standard deviation values to be applied for all grid cells, for each environmental tolerance threshold.
#' @param plasticity A list of length n.env.vars specifying the plasticity in the threshold tolerance trait for each environmental variable.
#' @return A set of output files written in the specified folder, each containing the run.name
#' @examples
#' ## Complete example of the multiple steps in applying AdaptR using data example provided for Drosophila jambulina
#'
#' # 1. Set the filepath to the example data set
#' filepath.data <- system.file("extdata", package="AdaptR")
#'
#' # 2. Create text files to describe the file path to each input variable, and the output variable
#' write.table((file.path(filepath.data,"Tmax",paste0("Tmax",seq(1:159),".asc"))), file = file.path(filepath.data,"Tmax","Tmax_filenames.txt"), eol = "\n", row.names = FALSE, col.names = FALSE, quote=FALSE )
#' write.table((file.path(filepath.data,"Habitat",paste0("Habitat",seq(1:159),".asc"))), file = file.path(filepath.data,"Habitat","Habitat_filenames.txt"), eol = "\n", row.names = FALSE, col.names = FALSE, quote=FALSE )
#' write.table((file.path(filepath.data,"compact_grids",paste0("demo_compact_grids_T",seq(1:159)))), file = file.path(filepath.data,"compact_grids","demo_compact_grids_output_filenames.txt"), eol = "\n", row.names = FALSE, col.names = FALSE, quote=FALSE )
#' 
#' # 3. Run the grid compactor
#' CompactGrids(compactor.parameter.file.name = file.path(filepath.data,"species_inputs","Jambulina__grid_compactor_parameter_file.txt"),
#'              ncols = 100,
#'              nrows = 79,
#'              n.env.vars = 2,
#'              n.time.points = 159,
#'              raw.env.grids.name.files = c(file.path(filepath.data,"Tmax","Tmax_filenames.txt"),file.path(filepath.data,"Habitat","Habitat_filenames.txt")),
#'              output.env.name.file = file.path(filepath.data,"compact_grids","demo_compact_grids_output_filenames.txt"))
#'
#' # 4. Generate a dispersal kernel for the drosophila example
#' filepath.data <- system.file("extdata", package="AdaptR")
#' Dispersal_Neighbourhood(radius=5, 
#'                        type="neg.power", 
#'                        params=c(1,1), 
#'                        output.name="Dispersal_relfile_L1_K1_rad5",
#'                        output.directory=file.path(filepath.data,"species_inputs"),
#'                        dispersal.plot=TRUE)  
#' 
#' # 5. Create text files to describe the file path to the compact grids
#' write.table(paste0("demo_compact_grids_T",rep(1:159, length.out=159)), file = file.path(filepath.data,"species_inputs","compact_series_names.txt"), eol = "\n", row.names = FALSE, col.names = FALSE, quote=FALSE )
#' 
#' # 6. Run the AdaptR model
#' AdaptR(run.name = "jambulina_test",
#'         parameter.file.name = file.path(filepath.data,"outputs","jambulina_test_parameters.txt"),
#'         ncols = 100,
#'         nrows = 79,
#'         output.folder.path = file.path(filepath.data,"outputs"),
#'         verbose.outputs = FALSE,
#'         n.time.points = 159,
#'         n.env.vars = 2,
#'         env.vars.names = c("MaxTemp", "Other_Maxent"),
#'         env.grids.folder.path = file.path(filepath.data,"compact_grids"),
#'         env.grids.name.file = file.path(filepath.data,"species_inputs","compact_series_names.txt"),
#'         species.initial.grid = file.path(filepath.data,"species_inputs","demo_jambulia_initial_distribution.asc"),
#'         minimum.survival.percentage = 5,
#'         resident.population.weighting = 1000,
#'         dispersal.neighbourhood.file = file.path(filepath.data,"species_inputs","Dispersal_relfile_L1_K1_rad5.dna"),
#'         species.location.file = file.path(filepath.data,"species_inputs","demo_locations_out.txt"),
#'         env.lower.thresholds = c(19.47,0.99),
#'         env.upper.thresholds = c(37.94547,1.01),
#'         env.low.adaptation = c(FALSE,FALSE),
#'         env.high.adaptation = c(TRUE,FALSE),
#'         adapt.limit = c(40,0),
#'         heritability = c(0.53,0),
#'         fitness.cost = c(0.05,0),
#'         adapt.threshold.grids = c(FALSE,FALSE),
#'         phenotypic.sd.grid = c(FALSE,FALSE),
#'         phenotypic.sd.value = c(1.106,0),
#'         plasticity = c(1.106,0))
#'
#' # 7. Map the predictions of the model (at the last time point)
#' map.single.run(output.folder.path = file.path(filepath.data,"outputs"), 
#'                run.name = "jambulina_test")
#'
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib AdaptR
#' @export
AdaptR <- 
function(run.name,
                 parameter.file.name,
                 ncols,
                 nrows,
                 output.folder.path,
                 verbose.outputs = FALSE,
                 n.time.points,
                 n.env.vars,
                 env.vars.names,
                 env.grids.folder.path,
                 env.grids.name.file,
                 species.initial.grid,
                 minimum.survival.percentage = 5,
                 resident.population.weighting = 1000,
                 dispersal.neighbourhood.file,
                 species.location.file,
                 env.lower.thresholds,
                 env.upper.thresholds,
                 env.low.adaptation,
                 env.high.adaptation,
                 adapt.limit,
                 heritability,
                 fitness.cost,
                 adapt.threshold.grids,
                 adapt.threshold.grid.names,
                 phenotypic.sd.grid,
                 phenotypic.sd.grid.names,
                 phenotypic.sd.value,
                 plasticity)
{
  ##_____________________________________________________________________________________##
  ## First, let's run some sanity checks on the parameters
  if(missing(run.name))
    stop("Need to specify a name for this run of Adaptor.")  
  if(!file.exists(output.folder.path))
    stop("Need to specify valid output file path.")
  output.folder.path<-paste0(output.folder.path,.Platform$file.sep)
  if(missing(ncols))
    stop("Need to specify the number of columns in the grid (ncols).")
  if(missing(nrows))
    stop("Need to specify the number of rows in the grid (nrows).")  
  if(missing(n.time.points))
    stop("Need to specify the number of time points to be modelled (n.time.points).")  
  if(missing(n.env.vars))
    stop("Need to specify the number of environmental variables in your compacted grids (n.env.vars).")  
  if(!(length(env.vars.names) == n.env.vars))
    stop("You have specified a different number of environment variable names to n.env.vars.")
  if(!file.exists(env.grids.folder.path))
    stop("Need to specify valid file path to the compacted environment grids.")
  env.grids.folder.path<-paste0(env.grids.folder.path,.Platform$file.sep)
  if(!file.exists(env.grids.name.file))
    stop("Need to specify valid file providing names for the environment grids.")  
  if(!file.exists(species.initial.grid))
    stop("Need to specify valid file for the species initial occurrence.") 
  if(!file.exists(dispersal.neighbourhood.file))
    stop("Need to specify valid file for the dispersal neighbuorhood probabillities.")  
  if(!file.exists(species.location.file))
    stop("Need to specify valid file for locations to track (species.location.file).")
  if(missing(env.lower.thresholds))
    stop("Need to specify lower thresholds for your environmental tolerance.") 
  if(length(env.lower.thresholds) < n.env.vars)
    stop(paste0("The lower tolerance thresholds need to be specified for all ",n.env.vars)," environmental variables.")   
  if(missing(env.upper.thresholds))
    stop("Need to specify upper thresholds for your environmental tolerance.") 
  if(length(env.upper.thresholds) < n.env.vars)
    stop(paste0("The upper tolerance thresholds need to be specified for all ",n.env.vars)," environmental variables.")   
  if(missing(env.low.adaptation))
    env.low.adaptation <- rep(FALSE,times=n.env.vars)   
  if(missing(env.high.adaptation))
    env.high.adaptation <- rep(FALSE,times=n.env.vars)  
  if(length(env.low.adaptation) < n.env.vars)
    env.low.adaptation <- c(env.low.adaptation,rep(FALSE,times=(n.env.vars-length(env.low.adaptation))))  
  if(length(env.high.adaptation) < n.env.vars)
    env.high.adaptation <- c(env.high.adaptation,rep(FALSE,times=(n.env.vars-length(env.high.adaptation))))       
  for(i.env in 1:n.env.vars)
  {
    if(env.low.adaptation[i.env] && env.high.adaptation[i.env] == TRUE)
      stop("You cannot specify adaptation on both the upper and lower tolerance threshold for the same variable.")
  }  
  if(missing(adapt.limit))
  {
    warning("No adaptation limits specified, so no adaptation will be implemented.")
    env.low.adaptation <- rep(FALSE,times=n.env.vars)
    env.high.adaptation <- rep(FALSE,times=n.env.vars)
  }
  if(!missing(adapt.limit))
  {
    if(length(adapt.limit) < n.env.vars)
    {
      warning("No adaptation limits are specified for some environment variables, so no adaptation will be implemented for those variables.")
      env.low.adaptation <- c(env.low.adaptation,rep(FALSE,times=(n.env.vars-length(adapt.limit))))
      env.high.adaptation <- c(env.high.adaptation,rep(FALSE,times=(n.env.vars-length(adapt.limit))))
    }
  }
  if(missing(heritability))
  {
    warning("No heritability values specified, so no adaptation will be implemented.")
    env.low.adaptation <- rep(FALSE,times=n.env.vars)
    env.high.adaptation <- rep(FALSE,times=n.env.vars)
  }
  if(!missing(heritability))
  {
    if(length(heritability) < n.env.vars)
    {
      warning("No heritability values are specified for some environment variables, so no adaptation will be implemented for those variables.")
      env.low.adaptation <- c(env.low.adaptation,rep(FALSE,times=(n.env.vars-length(heritability))))
      env.high.adaptation <- c(env.high.adaptation,rep(FALSE,times=(n.env.vars-length(heritability))))
    }
  }
  if(missing(fitness.cost))
  {
    warning("No fitness cost values specified, so none will be applied.")
    fitness.cost <- rep(0,times=n.env.vars)
  }
  if(length(fitness.cost)<n.env.vars)
  {
    warning("No fitness cost values specified for some variables, so none will be applied for them.")
    fitness.cost <- c(fitness.cost,rep(0,times=(n.env.vars-length(fitness.cost))))
  }    
  if(missing(adapt.threshold.grids))
  {
    adapt.threshold.grids <- rep(FALSE,times=n.env.vars)
  }
  if(length(adapt.threshold.grids)<n.env.vars)
  {
    adapt.threshold.grids <- c(dapt.threshold.grids,rep(FALSE,times=(n.env.vars-length(fitness.cost))))
  }    
  
  if(TRUE %in% adapt.threshold.grids)
  {
    if((missing(adapt.threshold.grid.names))||(length(adapt.threshold.grid.names)<n.env.vars))
    {
      stop("You need to specify valid filepaths to grids of initial threshold values under adaptation.") 
    }
    for(i.env in 1:n.env.vars)
    {
      if(adapt.threshold.grids[i.env] == TRUE)
      {
        if(!file.exists(adapt.threshold.grid.names[i.env]))
          stop("The filepath to your grids of initial threshold values under adaptation are not valid.")
      } #end if
    } # end for i.env     
  } # end if TRUE 
  if(TRUE %in% phenotypic.sd.grid)
  {
    if((missing(phenotypic.sd.grid.names))||(length(phenotypic.sd.grid.names)<n.env.vars))
    {
      stop("You need to specify valid filepaths to grids of initial phenotypic sd values under adaptation.") 
    }
    for(i.env in 1:n.env.vars)
    {
      if(phenotypic.sd.grid[i.env] == TRUE)
      {
        if(!file.exists(phenotypic.sd.grid.names[i.env]))
          stop("The filepath to your grids of initial phenotypic sd values under adaptation are not valid.")
      } #end if
    } # end for i.env     
  } # end if TRUE 
  
  if((TRUE %in% env.low.adaptation)||(TRUE %in% env.high.adaptation))
  {
    for(i.env in 1:n.env.vars)
    {
      if((env.low.adaptation[i.env] == TRUE) || (env.high.adaptation[i.env] == TRUE))
      { 
        if(phenotypic.sd.grid[i.env] == FALSE)
        {
          if(missing(phenotypic.sd.value))
            stop("You need to specify an initial value of phenotypic sd values under adaptation.")
          if(length(phenotypic.sd.value)<n.env.vars)  
            stop("You need to correctly specify an initial value of phenotypic sd values under adaptation.") 
        } # end if phenotypic.sd.grid.names[i.env] == FALSE
      } # end if
    } # end for i.env
  } # end if TRUE %in%
  
  
  ##_____________________________________________________________________________________##
  ## Now, write the parameter file...
  parameter.file = parameter.file.name                    ##CHANGED TO TAKE A SPECIFIED PARAMETER FILE NAME
  jobscript = file(parameter.file,'w')   # 'wb' might create a "unix-friendly" binary mode
  if(missing(plasticity))
  {
    plasticity <- rep(0,times=n.env.vars)
  }
  if(!missing(plasticity))
  {
    if(length(plasticity) < n.env.vars)
    {
      warning("No plasticity values are specified for some environment variables, so no plasticity will be implemented for those variables.")
      plasticity <- c(plasticity,rep(0,times=(n.env.vars-length(plasticity))))
    }
  }    
  writeLines(paste0("NCOLS ",ncols), con = jobscript)
  writeLines(paste0("NROWS ",nrows), con = jobscript)
  writeLines("OUTPUT_FOLDER", con = jobscript)
  writeLines(output.folder.path, con = jobscript)
  writeLines(paste0("OUTPUT_ROOTNAME ",run.name), con = jobscript)
  if(verbose.outputs)
  {
    writeLines("WRITE_EACH_STEP 1", con = jobscript)  
  }
  else
  {
    writeLines("WRITE_EACH_STEP 0", con = jobscript)
  }  
  writeLines(paste0("NUM_TIME_POINTS ",n.time.points), con = jobscript)
  writeLines(paste0("NUM_ENVIRONMENT_VARIABLES ", n.env.vars), con = jobscript)
  writeLines("ENV_GRIDS_FOLDER", con = jobscript)
  writeLines(env.grids.folder.path, con = jobscript)  
  writeLines("ENV_GRIDS_NAMEFILE", con = jobscript)
  writeLines(env.grids.name.file, con = jobscript)  
  writeLines("SPECIES_OCCURRENCE_GRID", con = jobscript)
  writeLines(species.initial.grid, con = jobscript)
  writeLines("MIN_SURVIVAL_PERCENTAGE", con = jobscript)
  writeLines(as.character(minimum.survival.percentage), con = jobscript)
  writeLines("RESIDENT_POPULATION_ATTRIBUTE_WEIGHTING", con = jobscript)
  writeLines(as.character(resident.population.weighting), con = jobscript)
  writeLines("DISPERSAL_NEIGHBOURHOOD_FILE", con = jobscript)
  writeLines(dispersal.neighbourhood.file, con = jobscript)
  writeLines("SPECIES_LOCATION_FILE", con = jobscript)
  writeLines(species.location.file, con = jobscript)
  # Now loop through each env variable, and write it's details to the parameter file
  for(i.env in 1:n.env.vars)
  {
    writeLines(paste0("##ENV_VARIABLE",i.env,"-",env.vars.names[i.env],"##"), con = jobscript)
    writeLines("LOWFUNDLIMIT", con = jobscript)
    writeLines(as.character(env.lower.thresholds[i.env]), con = jobscript)
    writeLines("HIGHFUNDLIMIT", con = jobscript)
    writeLines(as.character(env.upper.thresholds[i.env]), con = jobscript)
    writeLines("LOWADAPTATION", con = jobscript)
    if(env.low.adaptation[i.env])
    {
      writeLines("1", con = jobscript)
    }
    else
    {
      writeLines("0", con = jobscript)
    }
    writeLines("HIGHADAPTATION", con = jobscript)
    if(env.high.adaptation[i.env])
    {
      writeLines("1", con = jobscript)
    }
    else
    {
      writeLines("0", con = jobscript)
    }    
    if(env.low.adaptation[i.env] || env.high.adaptation[i.env])
    {
      writeLines("ADAPTATIONLIMIT", con = jobscript)
      writeLines(as.character(adapt.limit[i.env]), con = jobscript)
      writeLines("HERITABILITY", con = jobscript)
      writeLines(as.character(heritability[i.env]), con = jobscript)
      writeLines("FITNESSCOST", con = jobscript)
      writeLines(as.character(fitness.cost[i.env]), con = jobscript)
      if(adapt.threshold.grids[i.env]) 
      {
        writeLines("THRESHOLD_GRID", con = jobscript)
        writeLines("1", con = jobscript)
        writeLines("THRESHOLD_GRID_NAME", con = jobscript)
        writeLines(adapt.threshold.grid.names[i.env], con = jobscript)
      } # end if adapt.threshold.grids[i.env]
      else
      {
        writeLines("THRESHOLD_GRID", con = jobscript)
        writeLines("0", con = jobscript)
        writeLines("THRESHOLD_VALUE", con = jobscript)        
        if(env.low.adaptation[i.env])
        {
          writeLines(as.character(env.lower.thresholds[i.env]), con = jobscript)
        }
        else
        {
          writeLines(as.character(env.upper.thresholds[i.env]), con = jobscript)
        }
      } # end else adapt.threshold.grids[i.env]
      if(phenotypic.sd.grid[i.env])
      {
        writeLines("STD_DEV_GRID", con = jobscript)
        writeLines("1", con = jobscript)
        writeLines("PHENOTYPIC_SD_GRID_NAME", con = jobscript)
        writeLines(as.character(phenotypic.sd.grid.names[i.env]), con = jobscript)        
      }
      else
      {
        writeLines("STD_DEV_GRID", con = jobscript)
        writeLines("0", con = jobscript)
        writeLines("PHENOTYPIC_SD_VALUE", con = jobscript)
        writeLines(as.character(phenotypic.sd.value[i.env]), con = jobscript)
      } 
      writeLines("PLASTICITY", con = jobscript)
      writeLines(as.character(plasticity[i.env]), con = jobscript)
    } #end if env.low.adaptation[i.env] || env.high.adaptation[i.env]
  } #end for i.env
  writeLines("PERCENTAGE_SUCCESS 100",  con = jobscript)
  close(jobscript)
  
  ## Now test that the parameter file has been written out properly
  if(!file.exists(parameter.file))
    stop("Creation of the parameter file has failed. Ensure a correct file path and suitable name have been provided for the 'parameter.file.name' argument")  
  
  ##_____________________________________________________________________________________##
  # load the dll
#  package.path<-system.file(package="AdaptR")
#  r_arch <- .Platform$r_arch
#  file.path.source<-file.path(package.path, "libs", r_arch, "main.dll")
  # load the dll
#  dyn.load(file.path.source)
  # call the AdaptR function in the dll
#  AdaptR.out <- .C("AdaptR",  argv = as.character(c(parameter.file)), arg_i_catch = as.integer(c(-3,0)))
  # unload the dll
#  dyn.unload(file.path.source)
  ##_____________________________________________________________________________________##
  AdaptR.out <- rcpp_AdaptR(parameter.file)
  
  ##_____________________________________________________________________________________##    
  # Now provide some output
  if(AdaptR.out == -3 )
    stop("AdaptR has not run. Could not call the dll.")
  if(AdaptR.out == -2 )
    stop("AdaptR has not run because the parameter file is formatted incorrectly.")  
  if(AdaptR.out == -1 )
    stop("AdaptR has not run because the parameter file does not exist.")  
  
  message(paste0("AdaptR completed. The results are in the specified output folder."))
  
  # library(roxygen2)
  # library(devtools)
  # note: enter > document() 
  # to load the NAMESPACE documentation correctly
  # To install package from Github:
  #  library(devtools)
  #  install_github("KarelMokany/AdaptR")
  
} # end AdaptR function



