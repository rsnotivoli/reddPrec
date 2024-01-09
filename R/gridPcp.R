#' Gridded dataset creation
#' 
#' @description This function creates a gridded precipitation dataset from a station-based dataset.
#' @param prec matrix or data.frame containing the original (cleaned) precipitation data. Each column represents one station. The names of columns must coincide with the names of the stations.
#' @param grid SpatRaster. Collection of rasters representing each one of the predictors.
#' @param sts matrix or data.frame. A column "ID" (unique ID of stations) is required. The rest of the columns (all of them) will act as predictors of the model.
#' @param dates vector of class "Date" with all days of observations (yyyy-mm-dd).
#' @param crs character. Coordinates system in EPSG format (e.g.: "EPSG:4326").
#' @param coords vector of two character elements. Names of the fields in "sts" containing longitude and latitude.
#' @param coords_as_preds logical. If TRUE (default), "coords" are also taken as predictors.
#' @param neibs integer. Number of nearest neighbors to use.
#' @param thres numeric. Maximum radius (in km) where neighboring stations will be searched. NA value uses the whole spatial domain.
#' @param ncpu number of processor cores used to parallel computing.
#' @export
#' @importFrom doParallel registerDoParallel
#' @examples
#' \dontrun{
#' alt <- terra::rast(volcano, crs = 'EPSG:4326')
#' terra::ext(alt) <- c(-1,3,38,42)
#' lon <- terra::rast(cbind(terra::crds(alt),terra::crds(alt)[,1]),type='xyz',crs='EPSG:4326')
#' lat <- terra::rast(cbind(terra::crds(alt),terra::crds(alt)[,2]),type='xyz',crs='EPSG:4326')
#' dcoast <- terra::costDist(alt,target=min(terra::values(alt)))/1000
#' grid <- c(alt, lon, lat, dcoast)
#' names(grid) <- c('alt', 'lon', 'lat', 'dcoast')
#' 
#' prec <- round(matrix(rnorm(2*25, mean = 1.2, sd = 4), 2, 25), 1)+1
#' prec[prec<0] <- 0
#' colnames(prec) <- paste0('sts_',1:25)
#' sts <- data.frame(ID = paste0('sts_',1:25), as.data.frame(terra::spatSample(grid, 25)))
#' gridPcp(prec, grid, sts, 
#'         dates = seq.Date(as.Date('2023-04-01'),as.Date('2023-04-02'),by='day'),
#'         thres = NA, coords = c('lon','lat'),coords_as_preds = TRUE, 
#'         crs = 'EPSG:4326', neibs = 10, ncpu = 2)
#' 
#' r <- terra::rast(c('./pred/20230401.tif','./err/20230401.tif'))
#' }

gridPcp <- function (prec, grid, sts, dates, ncpu, thres, neibs,
                     coords, crs, coords_as_preds){
  
  # checks
  m <- match(colnames(prec), sts$ID)
  prec <- prec[,m]
  
  
  registerDoParallel(cores=ncpu)
  
  for(i in 1:nrow(prec)){
    message(paste0('[',Sys.time(),'] -', " Computing day ",dates[i]))
    predday(x = prec[i,], 
            grid = grid, 
            sts = sts[,-which(colnames(sts)=='ID')], 
            neibs = neibs, coords = coords, crs = crs,
            coords_as_preds = coords_as_preds,
            date = dates[i],
            thres = thres)
  }
  
  message(paste0('[',Sys.time(),'] -', " END"))
}
