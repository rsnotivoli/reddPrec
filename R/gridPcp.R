
#' This function creates gridded precipitation data

#' @param filled matrix or data.frame containing the complete precipitation data by stations. Each column represents one station. The names of columns have to be names of the stations.
#' @param sts matrix or data.frame containing the stations info. Must have at least four fields: ID: station identifier; ALT: altitude; X: Longitude in UTM projection (meters); and Y: Latitude in UTM projection (meters). Tabulation separated.
#' @param inidate object of class character in format 'YYYY-mm-dd' defining the first day of quality control process
#' @param enddate object of class character in format 'YYYY-mm-dd' defining the last day of quality control process
#' @param neibs integer number of nearest neighbours to use
#' @param thres maximum radius where neighbouring stations will be searched
#' @param intermediate when TRUE, one file per day will be written in subdirectory ./days with the information about intermediate estimates.
#' @param ncpu number of processor cores used to parallel computing. It uses all the available cores by default.


gridPcp <- function (filled, points, sts, inidate, enddate, ncpu, 
                     thres = NA, neibs = 10, intermediate = T){
  dir.create("./gridded/", showWarnings = F)
  
  inidate <- as.Date(inidate)
  enddate <- as.Date(enddate)
  datess <- seq.Date(inidate, enddate, by = "day")
  
  #checks
  m <- match(colnames(filled), sts$ID)
  if(length(which(is.na(m))) > 0) 
    stop('Stations ID do not coincide with dataset names')
  
  #creates a distance matrix with the names of sorted neighbours
  distanc <- t(Apply(as.matrix(points[,c('LAT','LON')]), 
                     margins = 1, 
                     AtomicFun = '.dist_near', 
                     y = sts[,c('LAT','LON','ID')], 
                     thres = thres,
                     ncores = ncpu)[[1]])
  rownames(distanc) <- points$ID
  
  #select neighbours based on neibs
  distanc <- distanc[,1:neibs]
  
  gridded <- t(Apply(list(as.matrix(filled), as.matrix(1:length(datess))),
                       margins = 1, AtomicFun = '.predday', 
                       distanc = distanc, points = points, sts = sts, 
                       datess = datess, neibs = neibs, 
                       intermediate = intermediate, ncores = ncpu)[[1]])
  colnames(gridded) <- as.character(points$ID)
  
  roundd <- function(x){as.integer(round(x, 1)*10)}
  gridded <- apply(gridded, 2, roundd)
  
  save(gridded, points, file = 'gridded_pcp.RData')
  
}
