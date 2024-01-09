#' Standardization
#'
#' Function to standardize predictions to observations
#'
#' @param obs vector of observations
#' @param sim vector of predictions
#' @param method selected method of standardization
#' @param window odd integer. Length of data considered for standardization
#' @noRd


standardization <- function(obs, sim, method = NULL, dates, window){
  
  days <- substr(dates,6,10)
  
  if(!is.null(window)) if(window %% 2 == 0) stop("Window must be an odd integer")
  
  pred_st <- rep(NA,length = length(sim))
  if(is.null(method)){
    return(pred_st)
  } else {
    if(!is.null(window)){
      for(i in 1:length(obs)){
        # select days within the window
        ini <- (i-(window-1)/2)
        if(ini<1) ini <- 1
        end <- (i+(window-1)/2)
        if(end>length(obs)) end <- length(obs)
        dd <- dates[ini:end]
        # get data from all days
        m <- which(days%in%substr(dd,6,10))
        wd <- which(dates[i]==dates[m])
        o <- obs[m]
        s <- sim[m]
        if(method == "quant") rr <- stand_qq(o, s)
        if(method == "ratio") rr <- stand_ratio(o, s)
        pred_st[i] <- rr[wd]
      }
    } else {
      if(method == "quant") pred_st <- stand_qq(obs, sim)
      if(method == "ratio") pred_st <- stand_ratio(obs, sim)
    }
  }
  return(pred_st)
}
