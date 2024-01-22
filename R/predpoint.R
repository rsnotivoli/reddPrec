#' Estimation of new daily precipitation values for a single day and gridpoint
#' 
#' @description This function uses the neighbouring observations to estimate new precipitation values in a single day.
#' @param can SpatVector. Candidate gridpoint
#' @param ref SpatVector. Observations
#' @param thres numeric. Distance threshold to find neighbors.
#' @param neibs number of nearest neighbours that will be used to estimate new values
#' @param covars formula. Names of predictors.
#' @importFrom terra distance
#' @importFrom stats as.formula glm binomial predict
#' @noRd
#' 

predpoint <- function(can,ref,thres,neibs,covars,n){
  #set nearest observations
  dd <- terra::distance(can,ref)/1000
  if(!is.na(thres)){ 
    dd <- dd[dd<thres]
  }
  ref <- ref[match(sort(dd)[1:neibs],dd)]
  if (max(ref$val) == 0) {
    pred <- err <- 0
  } else{
    # probability of ocurrence prediction
    rr <- as.data.frame(ref)
    rr$val[rr$val > 0] <- 1
    
    f <- as.formula(paste0('val ~ ',covars))
    fmtb <- suppressWarnings(
      glm(f,family = binomial(),data = rr)
    )
    
    pb <- round(predict(fmtb, newdata = as.data.frame(can), 
                        type = "response"),2)
    #amount prediction
    #rescaling
    rr <- as.data.frame(ref)
    MINc <- min(rr$val) - (as.numeric(quantile(rr$val, 0.50)) - as.numeric(quantile(rr$val, 0.25)))
    MINc <- ifelse(MINc<0,0,MINc)
    MAXc <- max(rr$val) + (as.numeric(quantile(rr$val, 0.75))-as.numeric(quantile(rr$val, 0.50)))
    RANGE <- as.numeric(MAXc - MINc)
    rr$val <- (rr$val - MINc) / RANGE
    
    fmt <- suppressWarnings(
      glm(f,family = quasibinomial(),data = rr)
    )
    p <- predict(fmt, newdata = as.data.frame(can),type = "response")
    p <- round((p * RANGE) + MINc, 2)
    
    # error calculation
    e <- sqrt(sum((rr$val - predict(fmt, type = 'response')) ^ 2)/(length(rr$val) - length(n)))
    e <- round((e * RANGE) + MINc, 2)
    
    #evaluating estimate
    if(pb <= 0.5) pred <- 0 else pred <- p
    err <- e
  }
  return(c(pred,err))
}
