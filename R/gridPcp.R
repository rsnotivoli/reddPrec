#' Gridded dataset creation
#' 
#' @description This function creates a gridded precipitation dataset from a station-based dataset.
#' @param prec matrix or data.frame containing the original (cleaned) precipitation data. Each column represents one station. The names of columns must coincide with the names of the stations.
#' @param grid SpatRaster. Collection of rasters representing each one of the predictors.
#' @param dyncovars SpatRasterDataset. Collection of variables acting as dynamic predictors (changing each day). Each dataset inside the sds object (corresponding to each variable) must have the same rasters as number of days to grid. See Details.
#' @param sts matrix or data.frame. A column "ID" (unique ID of stations) is required. The rest of the columns (all of them) will act as predictors of the model.
#' @param model_fun function. A function that integrates the statistical hybrid model (classification and regression). The default is learner_glm, which is the original model. Other models are also available (learner_rf and learner_xgboost). Users can create their functions with different models as well.
#' @param dates vector of class "Date" with all days of observations (yyyy-mm-dd).
#' @param crs character. Coordinates system in EPSG format (e.g.: "EPSG:4326").
#' @param coords vector of two character elements. Names of the fields in "sts" containing longitude and latitude.
#' @param coords_as_preds logical. If TRUE (default), "coords" are also taken as predictors.
#' @param neibs integer. Number of nearest neighbors to use.
#' @param thres numeric. Maximum radius (in km) where neighboring stations will be searched. NA value uses the whole spatial domain.
#' @param dir_name character. Name of the of the folder in which the data will be saved. Default NA uses the original names.
#' @param ncpu number of processor cores used to parallel computing.
#' @details
#' All the rasters provided in "grid" and "dyncovars" must have the same spatial characteristics (resolution, crs, extent, etc.).
#' The function estimates precipitation based on the nearest observations ("sts" and "prec") using as covariates the predictors contained in "grid" and "dyncovars". 
#' Predictors of "grid" are used for all days while those on "dyncovars" are selected depending on the day. 
#' For instance, to model the first day the algorithm considers all rasters in "grid" and only those corresponding to the first day of each variable.
#' @export
#' @importFrom doParallel registerDoParallel
#' @examples
#' \dontrun{
#' # fixed covariates (elevation, latitude, longitude)
#' alt <- terra::rast(volcano, crs = 'EPSG:4326')
#' terra::ext(alt) <- c(-1,3,38,42)
#' lon <- terra::rast(cbind(terra::crds(alt),terra::crds(alt)[,1]),type='xyz',crs='EPSG:4326')
#' lat <- terra::rast(cbind(terra::crds(alt),terra::crds(alt)[,2]),type='xyz',crs='EPSG:4326')
#' dcoast <- terra::costDist(alt,target=min(terra::values(alt)))/1000
#' grid <- c(alt, lon, lat, dcoast)
#' names(grid) <- c('alt', 'lon', 'lat', 'dcoast')
#' 
#' # Dynamic covariates (Variable 1, Variable 2, Variable 3)
#' foo <- alt
#' terra::values(foo) <- runif(length(terra::values(alt)))
#' dyncovars1 <- rep(foo, 7)
#' names(dyncovars1) <- paste('dynvar1.day',1:terra::nlyr(dyncovars1), 
#'                       sep = "_") # not use blank space!
#' dyncovars2 <- dyncovars1*0.0234
#' names(dyncovars2) <- paste('dynvar2.day',1:terra::nlyr(dyncovars2), 
#'                       sep = "_") # not use blank space!
#' dyncovars3 <- dyncovars1*10502
#' names(dyncovars3) <- paste('dynvar3.day',1:terra::nlyr(dyncovars3), 
#'                       sep = "_") # not use blank space!
#' dyncovars <- terra::sds(dyncovars1, dyncovars2, dyncovars3)
#'  
#' # precipitation and stations generation
#' set.seed(123)
#' prec <- round(matrix(rnorm(7*25, mean = 1.2, sd = 4), 7, 25), 1)+1
#' prec[prec<0] <- 0
#' colnames(prec) <- paste0('sts_',1:25)
#' sts <- data.frame(ID = paste0('sts_',1:25), as.data.frame(terra::spatSample(grid, 25)))
#' 
#' 
#' gridPcp(prec = prec, 
#'         grid = grid, 
#'         sts = sts, 
#'         model_fun = learner_glm, 
#'         dates = seq.Date(as.Date('2023-04-01'),as.Date('2023-04-07'),by='day'), 
#'         ncpu = 4, 
#'         thres = NA, 
#'         neibs = 15,
#'         coords = c('lon','lat'), 
#'         crs = 'EPSG:4326', 
#'         coords_as_preds = TRUE)
#' 
#' r <- terra::rast(c('./pred/20230401.tif','./err/20230401.tif'))
#' terra::plot(r)
#' }

gridPcp <- function (prec, grid, dyncovars = NULL, sts, model_fun, dates, ncpu, thres, neibs,
                     coords, crs, coords_as_preds, dir_name = NA){
  
  registerDoParallel(cores=ncpu)
  
  # single day
  if(length(dates) == 1){
    if(!is.numeric(prec)) stop("Input data must be numeric")
    predday(x = prec, 
            grid = grid, 
            sts = sts[,-which(colnames(sts)=='ID')],
            model_fun = model_fun,
            neibs = neibs, coords = coords, crs = crs,
            coords_as_preds = coords_as_preds,
            date = dates,
            thres = thres,
            dir_name = dir_name)
  }
  
  # multiple days
  if(length(dates) > 1){
    
    # checks
    m <- match(colnames(prec), sts$ID)
    prec <- prec[,m]
    
    for(i in 1:nrow(prec)){
      message(paste0('[',Sys.time(),'] -', " Computing day ",dates[i]))
      date <- dates[i]
      
      if(!is.null(dyncovars)){
        # add data from dynamic covariates ----
        if(any(length(dates) != terra::nlyr(dyncovars))) 
          stop("Length of dates must be the same of the number of dyncovars layers")
        terra::time(dyncovars) <- dates
        st <- terra::vect(sts, geom = coords, crs = crs)
        n <- length(terra::varnames(dyncovars))
        dcovs <- NA
        dgrid <- list()
        for(ij in 1:n){
          sel <- dyncovars[[ij]][[which(terra::time(dyncovars[[ij]]) %in% date)]]
          # saving rasters
          dgrid[[ij]] <- sel
          # extract values to stations
          dcovs <- cbind(dcovs, terra::extract(x = sel, y = st, ID = FALSE))
        }
        dcovs <- dcovs[, -1]
        dgrid <- do.call(c, dgrid)
        # merge values at stations / adding temporal spatial fields
        stsm <- cbind(sts, dcovs)
        gridm <- c(grid, dgrid)

      } else{
        stsm <- sts
        gridm <- grid
      }

      predday(x = prec[i,],
              grid = gridm,
              sts = stsm[, -which(colnames(stsm)=='ID')],
              model_fun = model_fun,
              neibs = neibs, coords = coords, crs = crs,
              coords_as_preds = coords_as_preds,
              date = date,
              thres = thres,
              dir_name = dir_name)
    }
  }
  
  
  message(paste0('[',Sys.time(),'] -', " END"))
}
