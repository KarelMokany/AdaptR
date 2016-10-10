#' Generate summary maps from AdpatR outputs
#'
#' @param output.folder.path The directory where the output files are stored
#' @param summary.name An output name for the sumary files created
#' @examples
#' # Using data example provided for Drosophila jambulina
#'
#' @importFrom raster raster
#' @export
map.summary <- 
  function(output.folder.path = getwd(), 
           summary.name = "summary")
  {
  ##_____________________________________________________________________________________##
  ## First, let's run some sanity checks on the parameters
  if(missing(output.folder.path))
    warning("No output folder path specified: using working directory.")
  if(!file.exists(output.folder.path))
    stop("Output folder path invalid.")
  ## Now grab the relevant filenames from the folder
output.folder.path<-"C:/Users/mok010/Karels Files/Mokany Files/My R applications/AdaptR/data/outputs" # FOR TESTING
  f.occ.filepaths <- list.files(paste0(output.folder.path),pattern="_final_occurrence.asc",full=T) 
    
    
    
    
    
    # check the files are there  
  final.occurrence.fname <- file.path(output.folder.path,paste0(run.name,"_final_occurrence.asc")) 
  if(!file.exists(final.occurrence.fname))
    stop("Output files with this run name were not found in the specified folder.")
  ## find adaptive env grids for this run...
  list.files(paste0(output.folder.path),pattern="_final_Z.asc",full=T) ## UP TO HERE - list all files with the run name
  list.files(paste0(output.folder.path),pattern= run.name|"_final_Z.asc",full=T)
  
  if(n.env.vars.adapt > 0)
    {
    for(i.env in 1:n.env.vars.adapt)
      {
      env.z.fname <- file.path(output.folder.path,paste0(run.name,"_Env_",i.env,"_final_Z.asc"))
      env.vp.fname <- file.path(output.folder.path,paste0(run.name,"_Env_",i.env,"_final_Vp.asc"))
      if(!file.exists(env.z.fname))
        stop("Output files (adaptive environments) with this run name were not found in the specified folder.")
      if(!file.exists(env.vp.fname))
        stop("Output files (adaptive environments) with this run name were not found in the specified folder.")
      } # end for i.env
    } # end if n.env.vars.adapt > 0   
  ##_____________________________________________________________________________________##
  
  final.occurrence.ras <- raster(final.occurrence.fname)
  plot(final.occurrence.ras,main="Occupancy",col=c("grey","black"), legend=F)

  if(n.env.vars.adapt > 0)
    {
    for(i.env in 1:n.env.vars.adapt)
      {
      env.z.fname <- file.path(output.folder.path,paste0(run.name,"_Env_",i.env,"_final_Z.asc"))
      env.vp.fname <- file.path(output.folder.path,paste0(run.name,"_Env_",i.env,"_final_Vp.asc"))
      env.z.ras <- raster(env.z.fname)
      env.vp.ras <- raster(env.vp.fname)
      plot(final.occurrence.ras, main=paste0("Mean tolerance for environment ",i.env), col=c("grey","black"),legend=F)
      plot(env.z.ras, col=colours,add=T)
      plot(final.occurrence.ras, main=paste0("Variance in tolerance for environment ",i.env), col=c("grey","black"),legend=F)
      plot(env.vp.ras, col=colours,add=T)
      } # end for i.env
    } # end if n.env.vars.adapt > 0
  
  ##_____________________________________________________________________________________##
  
  
    
  ################ ALEXs CODE BELOW ###########################################################################  
  ## This was written for many replicate runs - to average & plot the results #################################    
  # @importFrom raster raster
  # @importFrom doParallel ??notsure??
  #map.run.summary.func = function(outdir, filepath.data, cores=4)
  #{
  # Land area ID and template
  basemap = raster(paste0(filepath.data,"/Tmax/Tmax1.asc")) ; basemap[which(!is.na(basemap[1:length(basemap)]))]=1
  land = which(basemap[1:length(basemap)]==1)
  
  #-------------------------------------------------------------------------------------------------------------#
  # CYCLE OVER FINAL OCCURRENCE MAPS
  #-------------------------------------------------------------------------------------------------------------#
  
  iterations = length(list.files(paste0(outdir),pattern="_final_occurrence.asc"))
  
  if(iterations==1)
    {
    fdist   = raster(list.files(paste0(outdir),pattern="_final_occurrence.asc",full=T))
    fdist   = as.vector(t(as.matrix(fdist)))[land]
    avg.occ = basemap ; avg.occ[land] = fdist ; rm(fdist)
    #
    traitZ  = raster(list.files(paste0(outdir),pattern="_Env_1_final_Z.asc",full=T)) 
    traitZ  = as.vector(t(as.matrix(traitZ)))[land]
    avg.Z = basemap ; avg.Z[land] = traitZ ; rm(traitZ)
    #
    traitVP = raster(list.files(paste0(outdir),pattern="_Env_1_final_Vp.asc",full=T)) 
    traitVP = as.vector(t(as.matrix(traitVP)))[land]
    avg.VP  = basemap ; avg.VP[land] = traitVP ; rm(traitVP)
    } # end if iterations==1
  else 
    {
    if(iterations<=20)
      {
      print("Combining multiple iterations in loop.")
      # Loop
      for(i in 1:iterations)
        {
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        fdist  = raster(list.files(paste0(outdir),pattern="_final_occurrence.asc",full=T)[i]) 
        fdist   = as.vector(t(as.matrix(fdist)))[land]
        if(exists("Fr")){ Fr = cbind(Fr,fdist) } else { Fr = matrix(fdist,ncol=1) }
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        traitZ = raster(list.files(paste0(outdir),pattern="_Env_1_final_Z.asc",full=T)[i]) 
        traitZ = as.vector(t(as.matrix(traitZ)))[land]
        if(exists("Zr")){ Zr = cbind(Zr,traitZ) } else { Zr = matrix(traitZ,ncol=1) }
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        traitVP = raster(list.files(paste0(outdir),pattern="_Env_1_final_Vp.asc",full=T)[i]) 
        traitVP = as.vector(t(as.matrix(traitVP)))[land]
        if(exists("Vr")){ Vr = cbind(Vr,traitVP) } else { Vr = matrix(traitVP,ncol=1) }
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        } # end for i
      # Convert to map of average occurrence (in the final generation)
      FrSum = rowSums(Fr) ; avg.occ = basemap ; avg.occ[land] = FrSum ; rm(Fr)
      # Convert to map of average occurrence (in the final generation)
      ZrAvg = rowSums(Zr,na.rm=T) / FrSum ; avg.Z = basemap ; avg.Z[land] = ZrAvg ; ZrAvg[ZrAvg==0] = NA
      # Convert to map of average occurrence (in the final generation)
      Vr[is.na(Vr)] = 0 ; VrAvg = rowSums(Vr, na.rm=T) / FrSum ; avg.VP = basemap ; avg.VP[land] = VrAvg ; rm(Vr, VrAvg)
      } # end if iterations<=20
    else
      {
      # Parallel loop
      print(paste0("Combining multiple iterations in parallel (",cores," cores)."))
      library(doParallel)
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
      cl = makeCluster(cores)
      registerDoParallel(cl)
      Fr = foreach(run=1:Nrun, .combine=cbind, .inorder=T, .packages="raster") %dopar%  {
        fdist   = raster(paste("outputs/Australia/",spp,"/",evo,"/",spp,"_",timeline,"_",evo,"_",run,"_final_occurrence.asc",sep="")) ; fdist   = as.vector(t(as.matrix(fdist)))[land]  }
      stopCluster(cl)
      # Convert to map of average occurrence (in the final generation)
      FrSum = rowSums(Fr) ; avg.occ = basemap ; avg.occ[land] = FrSum ; rm(Fr)
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
      cl = makeCluster(cores)
      registerDoParallel(cl)
      Zr = foreach(run=1:Nrun, .combine=cbind, .inorder=T, .packages="raster") %dopar%  {
        traitZ = raster(paste("outputs/Australia/",spp,"/",evo,"/",spp,"_",timeline,"_",evo,"_",run,"_Env_1_final_Z.asc",sep="")) ; traitZ = as.vector(t(as.matrix(traitZ)))[land]  }
      stopCluster(cl)
      # Convert to map of average occurrence (in the final generation)
      ZrAvg = rowSums(Zr,na.rm=T) / FrSum ; avg.Z = basemap ; avg.Z[land] = ZrAvg ; ZrAvg[ZrAvg==0] = NA
      Zr[is.na(Zr)] = 0
      indx = max.col(Zr, ties.method='first') ; ZrMax = Zr[cbind(1:nrow(Zr), indx)]; ZrMax[ZrMax==0] = NA ; max.CTmax = basemap ; max.CTmax[land] = ZrMax ; rm(Zr, ZrAvg, ZrMax)
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
      cl = makeCluster(cores) ; registerDoParallel(cl)
      Vr = foreach(run=1:Nrun, .combine=cbind, .inorder=T, .packages="raster") %dopar%  {
        traitVP = raster(paste("outputs/Australia/",spp,"/",evo,"/",spp,"_",timeline,"_",evo,"_",run,"_Env_1_final_Vp.asc", sep = ""))
        traitVP = as.vector(t(as.matrix(traitVP)))[land]
      }
      stopCluster(cl)
      # Convert to map of average occurrence (in the final generation)
      Vr[is.na(Vr)] = 0 ; VrAvg = rowSums(Vr, na.rm=T) / FrSum ; avg.VP = basemap ; avg.VP[land] = VrAvg ; rm(Vr, VrAvg)
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
      }# end else iterations<=20
    } # end else iterations==1
  
  # Save as file
  #writeRaster(avg.occ, paste0(outdir,"/Mean_Final_Occupancy.asc", sep = ""),overwrite=T)
  #writeRaster(avg.Z,   paste0(outdir,"/Mean_Final_Z.asc", sep = ""),overwrite=T)
  #writeRaster(avg.VP,  paste0(outdir,"/Mean_Final_Vp.asc", sep = ""),overwrite=T)
  
  # In some cases no runs survive to the end and can't be mapped
  if(any(as.matrix(as.vector(avg.occ))>0))
    {
    
    #-------------------------------------------------------------------------------------------------------------#
    # WRITE PLOT TO FILE
    #-------------------------------------------------------------------------------------------------------------#
    
    colours = colorRampPalette(c("dodgerblue4","dodgerblue","deepskyblue","springgreen2","chartreuse","chartreuse3","yellowgreen","greenyellow","yellow","goldenrod1","goldenrod3","orange","orangered","orangered3","orangered4"))(900)
    
    #jpeg(paste(outdir,"/Mean_Final_Occupancy.jpg",sep=""), height=7900, width=10000, res=300)
    plot(basemap, main="Occupancy", col="grey85", legend=F)
    plot(avg.occ, col=colours,add=T)
    #dev.off()
    
    #jpeg(paste(outdir,"/Mean_Final_CTmax.jpg",sep=""), height=7900, width=10000, res=300)
    plot(basemap, main="Mean of Selected Trait", col="grey85", legend=F)
    plot(avg.Z, col=colours,add=T)
    #dev.off()
    
    #jpeg(paste(outdir,"/Mean_Final_Vp.jpg",sep=""), height=7900, width=10000, res=300)
    plot(basemap, main="Mean Variance of Selected Trait", col="grey85", legend=F)
    plot(avg.VP, col=colours,add=T)
    #dev.off()
    } # end if any(as.matrix(as.vector(avg.occ))>0)
  else 
    { print("No data in final generation to map.") }
  
  } # End map.summary function
