
#' This function performs a complete Quality Control over daily precipitation

#' @param prec matrix or data.frame containing the original precipitation data. Each column represents one station. The names of columns have to be names of the stations.
#' @param sts matrix or data.frame containing the stations info. Must have at least four fields: ID: station identifier; ALT: altitude; X: Longitude in UTM projection (meters); and Y: Latitude in UTM projection (meters). Tabulation separated.
#' @param inidate object of class character in format 'YYYY-mm-dd' defining the first day of quality control process
#' @param enddate object of class character in format 'YYYY-mm-dd' defining the last day of quality control process
#' @param neibs integer number of nearest neighbours to use
#' @param thres maximum radius where neighbouring stations will be searched
#' @param printmeta when TRUE, one file per day will be written in subdirectory ./meta.
#' @param ncpu number of processor cores used to parallel computing. It uses all the available cores by default.

qcPrec <- function (prec, sts, inidate, enddate, neibs = 10, thres = NA,
                    printmeta = TRUE, ncpu = availableCores()) 
{
  inidate <- as.Date(inidate)
  enddate <- as.Date(enddate)
  datess <- seq.Date(inidate, enddate, by = "day")
  
  #checks
  m <- match(colnames(prec), sts$ID)
  if(length(which(is.na(m))) > 0) 
    stop('Stations ID do not coincide with dataset names')
  
  #set datasets
  ori = prec
  ids = colnames(prec)
  #order stations like the columns in dataset
  sts <- sts[m,]
  
  #creates a distance matrix with the names of sorted neighbours
  distanc <- t(Apply(as.matrix(sts[,c('LAT','LON')]), 
                     margins = 1, 
                     AtomicFun = 'dist_near', 
                     y = sts[,c('LAT','LON','ID')], 
                     thres = thres,
                     ncores = ncpu)[[1]])
  rownames(distanc) <- sts$ID
  
  #First round of iterations
  prec = cbind(prec, 0) #adding the iter control at end
  it = 1 #iter count
  seguir = 1 #stop iter control
  while (seguir == 1) {
    print(paste0('[',Sys.time(),'] -', " Iteration ", it, " of quality control"))
    prec <- t(Apply(as.matrix(prec), margins = 1, AtomicFun = '.qcFirst', 
              distanc = distanc, it = it, sts = sts, neibs = 10,
              ncores = ncpu)[[1]])
    colnames(prec) <- c(as.character(sts$ID),'0')
    #increase iteration
    it = it + 1 
    #check if need more iterations
    if (sum(prec[, ncol(prec)]) == 0) {
      seguir = 0
    }
  }
  prec <- prec[,1:(ncol(prec)-1)]
  
  #last iteration
  print(paste0('[',Sys.time(),'] -', "Last iteration of quality control"))
  prec <- t(Apply(list(as.matrix(prec), as.matrix(ori), as.matrix(1:length(datess))),
                margins = 1, AtomicFun = '.qcLast', 
                distanc = distanc, it = it, sts = sts, datess = datess, 
                printmeta = printmeta, ncores = ncpu)[[1]])
  
  #save file
  save(prec, sts, file = "cleaned.RData")
}
