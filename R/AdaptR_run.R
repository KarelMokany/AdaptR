#' Run AdaptR
#'
#' @param Pfile A file path to the parameter file
#' @return A set of output files in the specified location
#' @examples
#' CompactGrids("C:/MyDocuments/AdaptR_ParameterFile.txt")
AdaptR<-function(run.name,
                 parameter.file.name,
                 ncols,
                 nrows,
                 output.folder.path = as.character(getwd()),
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
  ##if(!file.exists(out.folder))
  ##  stop("Need to specify valid output file path.")  
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
  ##if(!file.exists(env.grids.folder.path))
  ##  stop("Need to specify valid file path to the compacted environment grids.")
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
  dyn.load("src/main.dll")
  # call AdaptR dll
  AdaptR.out <- .C("AdaptR",  argv = as.character(c(parameter.file)), arg_i_catch = as.integer(c(-3,0)))
  # unload the dll
  dyn.unload("src/main.dll")
  ##_____________________________________________________________________________________##
  
  # Now provide some output
  if(AdaptR.out$arg_i_catch[1] == -3 )
    stop("AdaptR has not run. Could not call the dll.")
  if(AdaptR.out$arg_i_catch[1] == -2 )
    stop("AdaptR has not run because the parameter file is formatted incorrectly.")  
  if(AdaptR.out$arg_i_catch[1] == -1 )
    stop("AdaptR has not run because the parameter file does not exist.")  
  
  message(paste0("AdaptR simulation completed in ",arg_i_catch[1]," seconds. The results are in the specified output folder."))
  
} # end AdaptR function



