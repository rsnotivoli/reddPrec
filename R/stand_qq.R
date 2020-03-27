

stand_qq <- function(obs, pred){
  ww <- which((obs + pred) != 0)
  if(length(ww) < 5) stop("Not enough data to fit with QQ method")
  qm.fit <- fitQmap(obs[ww], pred[ww], method ='QUANT')
  w <- which(!is.na(pred))
  if(length(w) < length(pred)){
    xx <- doQmap(pred[w], qm.fit)
    pred[w] <- xx
  } else{
    pred <- doQmap(pred, qm.fit)
  }
  return(pred)
}
