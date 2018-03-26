
#' This function fills the gaps in original data series of daily precipitation

#' @param prec matrix or data.frame containing the original (cleaned) precipitation data. Each column represents one station. The names of columns have to be names of the stations.
#' @param sts matrix or data.frame containing the stations info. Must have at least four fields: ID: station identifier; ALT: altitude; X: Longitude in UTM projection (meters); and Y: Latitude in UTM projection (meters). Tabulation separated.
#' @param inidate object of class character in format 'YYYY-mm-dd' defining the first day of quality control process
#' @param enddate object of class character in format 'YYYY-mm-dd' defining the last day of quality control process
#' @param ncpu number of processor cores used to parallel computing. It uses all the available cores by default.
#' @param thres maximum radius where neighbouring stations will be searched
#' @param neibs integer number of nearest neighbours to use
#' @param intermediate when TRUE, one file per day will be written in subdirectory ./days with the information about intermediate estimates.
#' @param validation when TRUE, scores are computed between observed and predicted data and a set of text files and plots will be written in subdirectory ./val.

gapFilling <- function(prec, sts, inidate, enddate, ncpu = 8, thres = NA, neibs = 10,
                       intermediate = TRUE, validation = TRUE){

  m <- match(colnames(prec), sts$ID)
  if(length(which(is.na(m))) > 0) 
    stop('Stations ID do not coincide with dataset names')
  
  datess <- seq.Date(as.Date(inidate), as.Date(enddate), by = "day")

  #order stations like the columns in dataset
  sts <- sts[m,]
  
  print(paste0('[',Sys.time(),'] -', " Computing distances"))
  #creates a distance matrix with the names of sorted neighbours
  distanc <- t(Apply(as.matrix(sts[,c('LAT','LON')]), 
                     margins = 1, 
                     AtomicFun = '.dist_near', 
                     y = sts[,c('LAT','LON','ID')], 
                     thres = thres,
                     ncores = ncpu)[[1]])
  rownames(distanc) <- sts$ID
  
  print(paste0('[',Sys.time(),'] -', " Filling gaps"))
  pred <- t(Apply(list(as.matrix(prec), as.matrix(1:length(datess))),
                  margins = 1, AtomicFun = '.fillData', 
                  distanc = distanc, sts = sts, datess = datess, 
                  neibs = neibs, intermediate = intermediate, 
                  ncores = ncpu)[[1]])
  
  print(paste0('[',Sys.time(),'] -', " Standardizing final data series"))
  #monthly standardization
  pred3 <- pred
  colnames(pred3) <- colnames(prec)
  for (i in 1:ncol(prec)) {
    for(h in 1:12){
      w <- which(h == as.numeric(substr(datess, 6, 7)))
      ww <- which(prec[w,i] + pred[w,i] != 0)
      ww2 <- which(prec[w,i] + pred[w,i] == 0)
      pred3[w,i][ww2] <- 0
      pred3[w,i] <- pred[w,i] / 
        ((sum(pred[w,i][ww], na.rm = T) + 1) / 
           (sum(prec[w,i][ww], na.rm = T) + 1))
    }
  }
  
  if(intermediate){
    print(paste0('[',Sys.time(),'] -', " Writing files"))
    aa <- list.files('./days/', full.names = TRUE)
    for (i in 1:length(aa)) {
      d = read.table(aa[i], header = T, sep = "\t")
      d$final_pred <- round(pred3[i, ], 1)
      write.table(d, aa[i], quote = F, row.names = F, sep = "\t", na = "")
    }
  }
  
  
  if(validation){
  p3 <- pred3
  print(paste0('[',Sys.time(),'] -', " Computing scores"))
    for(i in 1:ncol(prec)){p3[which(is.na(prec[,i])),i] <- NA}
    elev <- seq(0,max(sts$ALT))
    .scores(obs = prec, sim = p3, alts = c(elev[seq(1, length(elev), length(elev)/10)],max(elev)), 
           dates = seq.Date(as.Date(inidate), as.Date(enddate), by = 'day'), est = sts)
  }
  
  #replace gaps by estimates
  for (i in 1:ncol(prec)) {
    w <- which(is.na(prec[, i]))
    if (length(w) > 0)
      prec[w, i] <- as.integer(round(pred3[w, i], 1))
  }
  filled <- prec
  save(filled, file = 'Filled.RData')
}

  
  