#' Estimating new values in original missing values data series of daily precipitation
#' 
#' @description This function uses the neighboring observations to estimate new precipitation values in those days and locations where no records exist.
#' @param prec matrix containing the original (cleaned) precipitation data. Each column represents one station. The names of columns must coincide with the names of the stations.
#' @param sts data.frame. A column "ID" (unique ID of stations) is required. The rest of the columns (all of them) will act as predictors of the model.
#' @param model_fun function. A function that integrates the statistical hybrid model (classification and regression). The default is learner_glm, which is the original model. Other models are also available (learner_rf and learner_xgboost). Users can create their functions with different models as well.
#' @param dates vector of class "Date" with all days of observations (yyyy-mm-dd).
#' @param crs character. Coordinates system in EPSG format (e.g.: "EPSG:4326").
#' @param coords vector of two character elements. Names of the fields in "sts" containing longitude and latitude.
#' @param coords_as_preds logical. If TRUE (default), "coords" are also taken as predictors.
#' @param neibs integer. Number of nearest neighbors to use.
#' @param thres numeric. Maximum radius (in km) where neighboring stations will be searched. NA value uses the whole spatial domain.
#' @param stmethod standardization method. 'quant' or 'ratio', see details.
#' @param window odd integer. Length of data considered for standardization
#' @param ncpu number of processor cores used to parallel computing.
#' @export
#' @importFrom foreach foreach %dopar%
#' @importFrom reshape2 melt
#' @importFrom reshape cast
#' @importFrom doParallel registerDoParallel
#' @details
#' After the gap filling, "stmethod" allows for an standardization of the predictions based on the observations.
#' It only works for daily data. For other timescales (monthly, annual) use "stmethod=NULL".
#' The "window" parameter is a daily-moving centered window from which data is collected for each year 
#' (i.e. a 15-day window on 16th January will take all predictions from 1st to 30th January of all years to standardize them
#' with their corresponding observations. Only standardized prediction of 16th January is returned. Process is repeated for all days).
#' 
#' @examples
#' \dontrun{
#' set.seed(123)
#' prec <- round(matrix(rnorm(30*50, mean = 1.2, sd = 6), 30, 50), 1)
#' prec[prec<0] <- 0
#' prec <- apply(prec, 2, FUN = function(x){x[sample(length(x),5)] <- NA; x})
#' colnames(prec) <- paste0('sts_',1:50)
#' sts <- data.frame(ID = paste0('sts_',1:50), lon = rnorm(50,0,1), 
#'                   lat = rnorm(50,40,1), dcoast = rnorm(50,200,50))
#'filled <- gapFilling(prec, sts, 
#'                     dates = seq.Date(as.Date('2023-04-01'),
#'                     as.Date('2023-04-30'),by='day'), 
#'                     stmethod = "ratio", thres = NA, coords = c('lon','lat'),
#'                     coords_as_preds = TRUE, crs = 'EPSG:4326', neibs = 10, 
#'                     window = 11, ncpu = 2)
#'str(filled)
#'summary(filled)
#'}
#' 

gapFilling <- function(prec, sts, model_fun = learner_glm, dates, stmethod = NULL, thres = NA, neibs = 10,
                       coords, crs, coords_as_preds = TRUE, window, ncpu = 2){

  if(is.null(colnames(prec))){
    message('Guessed dataset names in same order as Stations ID')
    colnames(prec) <- sts$ID
  }
  
  m <- match(colnames(prec), sts$ID)
  if(length(which(is.na(m))) > 0) 
    stop('Stations ID do not coincide with dataset names')
  
  # reorder stations
  sts <- sts[m,]
  
  message(paste0('[',Sys.time(),'] -', " Filling gaps"))

  registerDoParallel(cores=ncpu)
  
  j <- NULL
  a <- foreach(j = 1:nrow(prec), .combine=rbind, .export=c("fillData")) %dopar% {
    fillData(x = prec[j,], 
            sts = sts[,-which(colnames(sts)=='ID')],
            model_fun = model_fun,
            neibs = neibs,
            coords = coords,
            crs = crs,
            coords_as_preds = coords_as_preds,
            thres=thres
            )
  }
  a <- data.frame(date = sort(rep(dates, ncol(prec))), ID = rep(sts$ID, nrow(prec)), a)
  
  
  
  
  message(paste0('[',Sys.time(),'] -', " Standardizing final data series"))
  
  bb <- a[,c('date','ID','mod_pred')]
  prec_pred <- suppressMessages(reshape::cast(bb,date~ID)[,-1])
  prec_pred <- prec_pred[,match(colnames(prec),colnames(prec_pred))]
  
  pred <- foreach(j = 1:ncol(prec_pred), .combine=cbind, .export=c("standardization","stand_qq")) %dopar% {
    standardization(obs = prec[,j],
                    sim = prec_pred[,j],
                    method = stmethod, 
                    dates = dates,
                    window = window
                    )
  }
  colnames(pred) <- colnames(prec_pred)
  pred <- data.frame(date = dates, pred, check.names = FALSE)
  pred <- reshape2::melt(pred, id.vars = 'date')
  pred <- pred[order(pred$date),]
  
  a <- data.frame(a[,1:6], st_pred = pred$value, a[,7:8])
  
  message(paste0('[',Sys.time(),'] -', " END"))  
  
  return(a)
  
}

  
  
