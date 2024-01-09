#' Estimation of new daily precipitation values for a single day
#' 
#' @description This function uses the neighbouring observations to estimate new precipitation values in a single day.
#' @param x numeric vector with all the observations in the candidate day.
#' @param grid SpatRaster. Collection of rasters representing each one of the predictors.
#' @param sts matrix or data.frame. A column "ID" (unique ID of stations) is required. The rest of the columns (all of them) will act as predictors of the model.
#' @param neibs integer. Number of nearest neighbors to use.
#' @param crs character. Coordinates system in EPSG format (e.g.: "EPSG:4326").
#' @param coords vector of two character elements. Names of the fields in "sts" containing longitude and latitude.
#' @param coords_as_preds logical. If TRUE (default), "coords" are also taken as predictors.
#' @param date value of class "Date"
#' @param thres numeric. Maximum radius (in km) where neighboring stations will be searched. NA value uses the whole spatial domain.
#' @importFrom foreach foreach %dopar%
#' @importFrom terra vect rast crds writeRaster as.points
#' @noRd
#' 
predday <- function(x, grid, sts, neibs, coords, crs, coords_as_preds, date, thres) {
  
  dir.create('pred', showWarnings = FALSE)
  dir.create('err', showWarnings = FALSE)
  
  n <- names(sts)
  if(!coords_as_preds) n <- n[-match(coords,names(sts))]
  if(neibs<length(n)*2) 
    stop('The number of neighbors must be at least the double of predictors')
  
  n <- names(sts)
  ng <-names(grid)
  if(!any(n%in%ng)) stop("Names of predictors must be the same on 
                            grid and observations")
  
  if(!coords_as_preds) n <- n[-match(coords,names(sts))]
  covars <- paste(n, collapse='+') # predictors
  
  #check input data
  ww <- which(is.infinite(x) | is.na(x))
  if(length(ww) > 0){
    x <- x[-ww]
    sts <- sts[-ww,]
    }
  
  if(length(x) < neibs){ # if less data than neighbours, ends
    rp <- rast(cbind(crds(grid),NA), type="xyz", crs=crs)
    re <- rast(cbind(crds(grid),NA), type="xyz", crs=crs)
    names(rp) <- names(re) <- date
    writeRaster(rp, paste0('./pred/',gsub('-','',date),'.tif'),overwrite=TRUE)
    writeRaster(re, paste0('./err/',gsub('-','',date),'.tif'),overwrite=TRUE)
  }  
    
  if(max(x) == 0) { # if max data is 0, ends
    rp <- rast(cbind(crds(grid),0), type="xyz", crs=crs)
    re <- rast(cbind(crds(grid),0), type="xyz", crs=crs)
    names(rp) <- names(re) <- date
    writeRaster(rp, paste0('./pred/',gsub('-','',date),'.tif'),overwrite=TRUE)
    writeRaster(re, paste0('./err/',gsub('-','',date),'.tif'),overwrite=TRUE)
  } else{
    
    ref <- data.frame(sts,val = x)
    ref <- vect(ref, geom = coords, crs = crs, keepgeom = TRUE)
    grd <- as.points(grid)  
    
  # start evaluating data
  j <- NULL
  rr <- foreach(j = 1:length(grd), .combine=cbind) %dopar% {
    predpoint(can = grd[j], ref = ref, thres = thres, neibs = neibs,
              covars=covars)
  }
    rr <- t(rr)
    
    rp <- rast(cbind(crds(grid),rr[,1]), type="xyz", crs=crs)
    re <- rast(cbind(crds(grid),rr[,2]), type="xyz", crs=crs)
    names(rp) <- names(re) <- date
    writeRaster(rp, paste0('./pred/',gsub('-','',date),'.tif'),overwrite=TRUE)
    writeRaster(re, paste0('./err/',gsub('-','',date),'.tif'),overwrite=TRUE)
    
  }
  
}
