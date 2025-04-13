#' Q-Q standardization
#'
#' Function to standardize with a Quantile-Quantile approach
#'
#' @param o vector of observations
#' @param s vector of predictions
#' @importFrom qmap fitQmap doQmap
#' @noRd

stand_qq <- function(o, s){
  w0 <- which(s == 0)
  ww <- which((o + s) != 0)
  if(length(ww) < 5) return(s) else{
    qm.fit <- fitQmap(o[ww],
                      s[ww],
                      method = "QUANT",
                      wet.day = FALSE,
                      qstep = 0.1)
    # wet.day was not set before (creating fake automatic wet.day)
    # wet.day = 0 = FALSE same    
    w <- which(!is.na(s))
    if(length(w) < length(s)){
      xx <- doQmap(s[w], qm.fit)
      s[w] <- xx
    } else{
      s <- doQmap(s, qm.fit)
    }
    s[w0] <- 0
    return(round(s,2))
    }
}
