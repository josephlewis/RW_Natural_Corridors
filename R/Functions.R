create_nc_surface <- function(cs, cells, dem, filename, ncores = 1) { 
  
  RcppParallel::setThreadOptions(numThreads = ncores)
  
  conductanceMatrix <- summary(cs$conductanceMatrix)
  conductanceMatrix$x <- 1 / conductanceMatrix$x
  directed_graph <- cppRouting::makegraph(conductanceMatrix,directed=TRUE)
  rm(cs)
  rm(conductanceMatrix)
  gc()
  
  message("created directed graph")
  
  dem_empty <- dem
  dem_empty[!is.na(dem_empty)] <- 0
  
  df_list <- split(cells, ceiling(seq_along(cells)/ncores))
  
  for(i in 1:length(df_list)) {
    
    message(paste0(i, " of ", length(df_list)))

    lcps <- cppRouting::get_multi_paths(Graph=directed_graph, from = df_list[[i]], to = cells, long = TRUE)
    
    lcps_table <- tabulate(as.numeric(lcps[,3]))
    lcps_nodes <- which(lcps_table != 0)
    lcps_vals <- lcps_table[lcps_nodes]
    
    dem_empty[lcps_nodes] <- dem_empty[lcps_nodes] + lcps_vals
    
  }
  
  terra::writeRaster(dem_empty, paste0("./Output/Natural_Corridors/", filename , ".tif"))
  
}

normalise_raster <- function(x){(x-terra::minmax(x)[1])/(terra::minmax(x)[2]-terra::minmax(x)[1])}