#' Generate summary maps from a single run of AdpatR
#'
#' @param output.folder.path The directory where the output files are stored
#' @param run.name The name for the AdaptR simulation run
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
#' @importFrom raster raster
#' @export
map.single.run <- 
  function(output.folder.path = getwd(), 
           run.name)
  {
  ##_____________________________________________________________________________________##
  ## First, let's run some sanity checks on the parameters
  if(missing(run.name))
    stop("Need to specify a name for this run of Adaptor.")  
  if(missing(output.folder.path))
    warning("No output folder path specified: using working directory.")
  if(!file.exists(output.folder.path))
    stop("Output folder path invalid.")  # check the files are there  
  final.occurrence.fname <- file.path(output.folder.path,paste0(run.name,"_final_occurrence.asc")) 
  if(!file.exists(final.occurrence.fname))
    stop("Output files with this run name were not found in the specified folder.")
  ## find adaptive env grids for this run...
  Z.files<-list.files(paste0(output.folder.path), pattern=glob2rx(paste0(run.name,"*_final_Z.asc")),full=T)
  Vp.files<-list.files(paste0(output.folder.path), pattern=glob2rx(paste0(run.name,"*_final_Vp.asc")),full=T)
  if(length(Z.files) != length(Vp.files))
    stop("Different number of output grids for env tolerance mean (z) and stdev (vp).")
  
  if(length(Z.files) > 0)
    {
    for(i.env in 1:length(Z.files))
      {
      env.z.fname <- Z.files[i.env]
      env.vp.fname <- Vp.files[i.env]
      if(!file.exists(env.z.fname))
        stop("Output files (adaptive environments) with this run name were not found in the specified folder.")
      if(!file.exists(env.vp.fname))
        stop("Output files (adaptive environments) with this run name were not found in the specified folder.")
      } # end for i.env
    } # end if n.env.vars.adapt > 0   
  ##_____________________________________________________________________________________##
  
  final.occurrence.ras <- raster(final.occurrence.fname)
  plot(final.occurrence.ras,main="Occupancy",col=c("grey","black"), legend=F)

  if(ength(Z.files) > 0)
    {
    for(i.env in 1:length(Z.files))
      {
      env.z.fname <- Z.files[i.env]
      env.vp.fname <- Vp.files[i.env]
      env.z.ras <- raster(env.z.fname)
      env.vp.ras <- raster(env.vp.fname)
      env.ID <- substr(Z.files[i.env], nchar(Z.files[i.env])-13+1, nchar(Z.files[i.env])-12)
      plot(final.occurrence.ras, main=paste0("Mean tolerance for environment ",env.ID), col=c("grey","black"),legend=F)
      plot(env.z.ras, col=colours,add=T)
      plot(final.occurrence.ras, main=paste0("Variance in tolerance for environment ",env.ID), col=c("grey","black"),legend=F)
      plot(env.vp.ras, col=colours,add=T)
      } # end for i.env
    } # end if n.env.vars.adapt > 0
  
  ##_____________________________________________________________________________________##

  } # End map.single.run function
