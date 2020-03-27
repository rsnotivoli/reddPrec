
#' This function fills the gaps in original data series of daily precipitation

#' @param prec matrix or data.frame containing the original (cleaned) precipitation data. Each column represents one station. The names of columns have to be names of the stations.
#' @param sts matrix or data.frame containing the stations info. Must have at least four fields: ID: station identifier; ALT: altitude; X: Longitude in UTM projection (meters); and Y: Latitude in UTM projection (meters). Tabulation separated.
#' @param inidate object of class character in format 'YYYY-mm-dd' defining the first day of quality control process
#' @param enddate object of class character in format 'YYYY-mm-dd' defining the last day of quality control process
#' @param stmethod standardization method. 'quant' or 'ratio', see details.
#' @param ncpu number of processor cores used to parallel computing. It uses all the available cores by default.
#' @param thres maximum radius where neighbouring stations will be searched
#' @param neibs number of nearest neighbours that will be used to estimate new values
#' @param intermediate when TRUE, one file per day will be written in subdirectory ./days with the information about intermediate estimates.
#' @param validation when TRUE, scores are computed between observed and predicted data and a set of text files and plots will be written in subdirectory ./val.

gapFilling <- function(prec, sts, inidate, enddate, stmethod = NULL, ncpu = 8, thres = NA, neibs = 10,
                       intermediate = TRUE, validation = TRUE){

  if(is.null(colnames(prec))){
    message('Guessed dataset names in same order as Stations ID')
    colnames(prec) <- sts$ID
  }
  
  m <- match(colnames(prec), sts$ID)
  if(length(which(is.na(m))) > 0) 
    stop('Stations ID do not coincide with dataset names')
  
  datess <- seq.Date(as.Date(inidate), as.Date(enddate), by = "day")

  # reorder stations
  sts <- sts[m,]
  
  message(paste0('[',Sys.time(),'] -', " Computing distances"))
  # creates a distance matrix with the names of sorted neighbours
  nlim <- nrow(sts)
  distanc <- t(Apply(as.matrix(sts[,c('LAT','LON')]), 
                     margins = 1, 
                     fun = 'dist_near', 
                     y = sts[,c('LAT','LON','ID')], 
                     thres = thres,
                     nlim = nlim,
                     ncores = ncpu)[[1]])
  rownames(distanc) <- sts$ID
  
  message(paste0('[',Sys.time(),'] -', " Filling gaps"))
  pred <- t(Apply(list(as.matrix(prec), as.matrix(1:length(datess))),
                  margins = 1, fun = 'fillData', 
                  distanc = distanc, sts = sts, datess = datess, 
                  neibs = neibs, intermediate = intermediate, 
                  ncores = ncpu)[[1]])
  
  message(paste0('[',Sys.time(),'] -', " Standardizing final data series"))
  
  # standardization
  pred3 <- Apply(list(as.matrix(prec), as.matrix(pred)), 
        margins = 2,
        fun = 'standardization',
        method = stmethod,
        dates = datess,
        ncores = ncpu)[[1]]
  colnames(pred3) <- colnames(prec)
  
  if(intermediate){
    message(paste0('[',Sys.time(),'] -', " Writing files"))
    aa <- list.files('./days/', full.names = TRUE)
    for (i in 1:length(aa)) {
      d = read.table(aa[i], header = T, sep = "\t")
      d$std_pred <- round(pred3[i, ], 1)
      write.table(d, aa[i], quote = F, row.names = F, sep = "\t", na = "")
    }
  }
  
  
  #replace gaps by estimates
  for (i in 1:ncol(prec)) {
    w <- which(is.na(prec[, i]))
    if (length(w) > 0)
      prec[w, i] <- as.integer(round(pred3[w, i], 1))
  }
  filled <- prec
  save(filled, file = 'Filled.RData')
  
  
  if(validation){
  p3 <- pred3
  message(paste0('[',Sys.time(),'] -', " Computing scores"))
    for(i in 1:ncol(prec)) p3[which(is.na(prec[,i])),i] <- NA
    elev <- seq(0,max(sts$ALT))
    scores(obs = prec, 
           sim = p3, 
           alts = c(elev[seq(1, length(elev), length(elev)/10)],max(elev)), 
           dates = seq.Date(as.Date(inidate), as.Date(enddate), by = 'day'), 
           est = sts)
  }
}

  
  
