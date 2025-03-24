#' Estimation of new daily precipitation values for a single day and gridpoint
#' 
#' @description This function uses the neighbouring observations to estimate new precipitation values in a single day.
#' @param can SpatVector. Candidate gridpoint
#' @param ref SpatVector. Observations
#' @param model_fun function. A function that integrates the statistical hybrid model (classification and regression)
#' @param thres numeric. Distance threshold to find neighbors.
#' @param neibs number of nearest neighbours that will be used to estimate new values
#' @param covars formula. Names of predictors.
#' @importFrom terra distance
#' @importFrom stats as.formula glm binomial predict
#' @noRd
#' 

predpoint <- function(can, ref, model_fun, thres, neibs, covars){
  #set nearest observations
  dd <- terra::distance(can,ref)/1000
  if(!is.na(thres)){ 
    dd <- dd[dd<thres]
  }
  ref <- ref[match(sort(dd)[1:neibs],dd)]
  if (max(ref$val) == 0) {
    pred <- err <- 0
  } else if (sum(diff(ref$val))==0){
    pred <- ref$val[1]
    err <- 0
    } else{
      
      out <- model_fun(ref = ref, can = can, covars = covars)
      out <- round(out, 2)
      pb <- out[1]
      p <- out[2]
      e <- out[3]
    
    #evaluating estimate
    if(pb <= 0.5) pred <- 0 else pred <- p
    err <- e
  }
  return(c(pred, err))
}
