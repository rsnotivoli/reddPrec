
#' Quality Control of daily precipitation observations

#' @description This function apply several threshold-based criteria to filter original observations of daily precipitation.
#' @param prec matrix or data.frame containing the original precipitation data. Each column represents one station. The names of columns have to be names of the stations.
#' @param sts matrix or data.frame. A column "ID" (unique ID of stations) is required. The rest of the columns (all of them) will act as predictors of the model.
#' @param crs character. Coordinates system in EPSG format (e.g.: "EPSG:4326").
#' @param coords vector of two character elements. Names of the fields in "sts" containing longitude and latitude.
#' @param coords_as_preds logical. If TRUE (default), "coords" are also taken as predictors.
#' @param neibs integer. Number of nearest neighbors to use.
#' @param thres numeric. Maximum radius (in km) where neighboring stations will be searched. NA value uses the whole spatial domain.
#' @param qc vector of strings with the QC criteria to apply. Default is "all". See details.
#' @param qc3 numeric. Indicates the threshold (number of times higher or lower) from which a observation, in comparison with its estimate, should be deleted. Default is 10.
#' @param qc4 numeric vector of length 2. Thresholds of wet probability (0 to 1) and magnitude (in the units of input precipitation data) from which a observation of value zero, in comparison with its estimate, should be deleted. Default is c(0.99, 5). 
#' @param qc5 numeric vector of length 2. Thresholds of dry probability (0 to 1) and magnitude (in the units of input precipitation data) from which a observation higher than a specific value (also in the original units), in comparison with its estimate, should be deleted. Default is c(0.01, 0.1, 5). 
#' @param ncpu number of processor cores used to parallel computing.
#' @export
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach
#' @details 
#' Parameter "sts" must have an "ID" field containing unique identifiers of the stations.
#' 
#' "qc" can be "all" (all criteria are applied) or a vector of strings (e.g.: c("1","2","4")) indicating the QC criteria to apply to observations:
#' "1" (suspect value): obs==0 & all(neibs>0); 
#' "2" (suspect zero): obs>0 & all(neibs==0); 
#' "3" (suspect outlier): obs is "qc3" times higher or lower than the estimate; 
#' "4" (suspect wet): obs==0 & wet probability > "qc4\[1\]" & estimate > "qc4\[2\]"; 
#' "5" (suspect dry): obs>"qc5\[3\]" & dry probability < "qc5\[1\]" & estimate < "qc5\[2\]"
#' @examples
#' \dontrun{
#' set.seed(123)
#' prec <- round(matrix(rnorm(30*50, mean = 1.2, sd = 6), 30, 50), 1)
#' prec[prec<0] <- 0
#' colnames(prec) <- paste0('sts_',1:50) 
#' sts <- data.frame(ID = paste0('sts_',1:50), lon = rnorm(50,0,1), 
#'                   lat = rnorm(50,40,1), dcoast = rnorm(50,200,50))
#' qcdata <- qcPrec(prec, sts, crs = 'EPSG:4326', coords = c('lon','lat'),
#'                  coords_as_preds = TRUE, neibs = 10, thres = NA,
#'                  qc = 'all', qc3 = 10, qc4 = c(0.99, 5), qc5 = c(0.01, 0.1, 5),
#'                  ncpu=2)
#'str(qcdata)
#'}

qcPrec <- function (prec, sts, crs, coords, coords_as_preds = TRUE, neibs = 10, thres = NA,
                    qc = 'all', qc3 = 10, qc4 = c(0.99, 5), qc5 = c(0.01, 0.1, 5), ncpu = 1) 
{
  
  #checks
  if(ncol(prec) != nrow(sts)) 
    stop('The number of stations and precipitation data is different')
  
  # same order of columns and stations
  prec <- prec[,match(sts$ID,colnames(prec))]
  
  registerDoParallel(cores=ncpu)
  
  #First round of iterations
  a <- cbind(prec, 0) #adding the iter control at end
  j <- NULL
  it <- 1 #iter count
  seguir <- 1 #stop iter control
  while (seguir == 1) {
    message(paste0('[',Sys.time(),'] -', " Iteration ", it, " of quality control"))
    
    a <- foreach(j = 1:nrow(a), .combine=cbind, .export=c("qcFirst")) %dopar% {
         qcFirst(x = a[j,], 
              it = it, 
              sts = sts[,-which(colnames(sts)=='ID')], 
              neibs = neibs,
              coords = coords,
              crs = crs,
              coords_as_preds = TRUE,
              thres = thres,
              qc = qc, qc3 = qc3, qc4 = qc4, qc5 = qc5)
    }
    a <- t(a)
    
    
    #increase iteration
    it <- it + 1 
    #check if need more iterations
    if (sum(a[, ncol(a)]) == 0) {
      seguir <- 0
    }
  }
  a <- a[,1:(ncol(a)-1)]
  rownames(a) <- NULL
  
  #last iteration
  message(paste0('[',Sys.time(),'] -', "Last iteration of quality control"))
  
  b <- foreach(j = 1:nrow(a), .combine=cbind, .export=c("qcLast")) %dopar% {
    qcLast(x = a[j,], 
           y = prec[j,],
           sts = sts[,-which(colnames(sts)=='ID')], 
           neibs = neibs,
           coords = coords,
           crs = crs,
           coords_as_preds = TRUE,
           thres = thres,
           qc = qc, qc3 = qc3, qc4 = qc4, qc5 = qc5)
  }
  cleaned <- t(b[1:nrow(sts),])
  rownames(cleaned) <- NULL
  codes <- t(b[(nrow(sts)+1):nrow(b),])
  rownames(codes) <- NULL
  
  message(paste0('[',Sys.time(),'] -', " End"))
  
  return(list(cleaned=cleaned, codes=codes))
}
