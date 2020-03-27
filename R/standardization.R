


standardization <- function(obs, pred, method = NULL, dates){
  pred_st <- rep(NA,length = length(pred))
  if(is.null(method)){
    return(pred_st)
  } else if(method == "quant"){
    for(i in 1:12){
      w <- which(i == as.numeric(substr(dates, 6, 7)))
      # If error, returns predicted values without standardization
      pred_st[w] <- tryCatch(stand_qq(obs[w], pred[w]),
      error = function(e) {return(pred[w])})
    }
  } else if(method == "ratio"){
    for(i in 1:12){
      w <- which(i == as.numeric(substr(dates, 6, 7)))
      # If error, returns predicted values without standardization
      pred_st[w] <- tryCatch(stand_ratio(obs[w], pred[w]),
      error = function(e) {return(pred[w])}      )
    }
  }
  return(pred_st)
}
