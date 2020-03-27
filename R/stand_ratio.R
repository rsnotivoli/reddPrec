
stand_ratio <- function(obs, pred){
  ww <- which((obs + pred) != 0)
  if(length(ww) < 5) stop("Not enough data to fit with RATIO method")
  pred <- pred / ((sum(pred[ww], na.rm = T) + 1) / (sum(obs[ww], na.rm = T) + 1))
  return(pred)
}
  