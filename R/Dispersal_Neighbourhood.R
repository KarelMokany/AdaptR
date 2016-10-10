#' Generate a dispersal neighbourhood file for AdaptR
#'
#' @param radius the maximum radius of the dispersal neighbourhood, in grid cell units
#' @param type the type of dispersal probability function to apply ('neg.power','neg.exponential')
#' @param params a vector of parameters for the dispersal function
#' @param output.name name for the dispersal neighbourhood 
#' @param output.directory output folder in which to save the dispersal neighbourhood file
#' @param dispersal.plot whether to display a plot of the resultant dispersal probability against distance
#' @return A dispersal neighbourhood for AdaptR output file (.dna) written to the specified folder, with the file labelled by the output.name
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
#' @export
Dispersal_Neighbourhood <- 
function(radius,
         type="neg.power",
         params=c(1,1),
         output.name,
         output.directory,
         dispersal.plot="FALSE")
  {
  ##_____________________________________________________________________________________##
  ## First, let's run some sanity checks on the parameters
  ## First, let's run some sanity checks on the parameters
  # test radius exists
  if(missing(radius))
    stop("Need to specify the dispesal radius, in grid cell units (radius).")  
  # test radius >1
  if(radius<1)
    stop("Dispersal radius needs to be '1' or greater.")  
  # test output name
  if(missing(output.name))
    stop("Need to specify a name for the dispersal neighbourhood (output.name).")  
  # test output directory
  if(missing(output.directory))
    stop("Need to specify a valid filepath for the dispersal neighbourhood file (output.directory).")  
  if(!file.exists(output.directory))
    stop("Need to specify valid filepath for the dispersal neighbourhood file (output.directory).")  

  ##_____________________________________________________________________________________##
  # set the parameters
  lambda<-params[1]
  if(length(params)>1)
    {K<-params[2]}
  
  # just need to loop through all the cells in the square neighbourhood
  disp.data <- matrix(NA, nrow=((((radius*2)+1))^2), ncol=3)
  start <- (-1*radius)
  end <- (radius)  
  n.rad.cells <- 0
  for(i.row in start:end)
    {
    for(i.col in start:end)
      {
      cell.dist<-sqrt((i.row^2)+(i.col^2))
      if(cell.dist <= radius)
        {
        if(type == "neg.power")
          {disp.prob <- (2*K)/((2*pi*cell.dist)*pi*lambda*(1+(cell.dist/lambda)*(cell.dist/lambda)))}
        if(type == "neg.exponential")
          {disp.prob <- (1/(2*pi*(lambda^2)))*exp((-1*cell.dist)/lambda)}
        n.rad.cells <- n.rad.cells + 1
        disp.data[n.rad.cells,1] <- i.row
        disp.data[n.rad.cells,2] <- i.col
        disp.data[n.rad.cells,3] <- disp.prob
        } # end if cell.dist
      }# end for i.col
    } #end for i.row
  
  # Clear out the empty part of the matrix & catch its length
  disp.data <- disp.data[complete.cases(disp.data),]
  n.cells <- nrow(disp.data)
  
  # round down the dispersal probs to '1' if they are greater than this
  disp.data[(disp.data[,3]>1),3] <- 1
  
  # now write the dispersal data to file
  dispersal.file <- file.path(output.directory,paste0(output.name,".dna"))
  jobscript = file(dispersal.file,'w')  
  writeLines(as.character(as.integer(n.cells)), con = jobscript)
  for(i.row in 1:nrow(disp.data))
    {
    writeLines(paste0(as.character(as.integer(disp.data[i.row,1]))," ",as.character(as.integer(disp.data[i.row,2]))," ",as.character(disp.data[i.row,3])), con = jobscript)
    }# end for i.row
  close(jobscript)
  
  ## Now test that the dispersal file has been written out properly
  if(!file.exists(dispersal.file))
    stop("Creation of the dispersal file has failed. Ensure a correct file path and suitable name have been provided for the 'output.name' and 'output.directory' arguments, so a dispersal file can be created.")  
  
  # Plot if specified  
  if(dispersal.plot)
  {
  plot(1, type="n", xlim=c(0,radius), ylim=c(0,max(disp.data[,3])), xlab = "distance from source cell (grid cells)", ylab="Probability of successful dispersal")
  x.dat<-seq(0,radius,length.out = 1000)
  if(type == "neg.power")
  {y.dat <- (2*K)/((2*pi*x.dat)*pi*lambda*(1+(x.dat/lambda)^2))}
  if(type == "neg.exponential")
    {y.dat <- (1/(2*pi*(lambda^2)))*exp((-1*x.dat)/lambda)}  
  y.dat[(y.dat>1)] <- 1
  lines(x.dat,y.dat,col="red") 
  } # end if dispersal.plot
  
  message(paste0("Dispersal Neighbourhood for AdaptR created. The .dna file created is in the specified output folder"))
  
  } # end function Dispersal_Neighbourhood
           