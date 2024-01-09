#' Ratio standardization
#'
#' Function to standardize with a ratio approach
#'
#' @param o vector of observations
#' @param s vector of predictions
#' @noRd
#' 
stand_ratio <- function(o, s){
  w0 <- which(s == 0)
  ww <- which((o + s) != 0)
  if(length(ww) < 5) return(s) else{
    s <- s / ((sum(s[ww], na.rm = T) + 0.1) / (sum(o[ww], na.rm = T) + 0.1))
    s[w0] <- 0
    return(round(s,2))
  }
}
  