hmg_cleaning_trend_ar1 <- function(time_serie, p_value = 0.05){
  
  
  # Apply detrending
  sample_value <- as.numeric(time_serie)
  mk_res <- Kendall::MannKendall(time_serie)
  
  # removing trend if mk_res$sl < p_value
  if (as.numeric(mk_res$sl) < p_value) {
    
    res <- pracma::detrend(sample_value, tt = "linear")
    zoo::coredata(time_serie) <- res
    
  } 
  
  #  Apply AR(1) removal
  sample_value <- as.numeric(time_serie)
  dwt_res <- car::durbinWatsonTest(stats::lm(sample_value ~ 1), max.lag = 1)
  
  # removing AR1 if dwt_res$p < p_value
  if (dwt_res$p < p_value) {
    
    ar_model <- stats::arima(sample_value, order = c(1, 0, 0))
    res <- as.numeric(ar_model$residual)
    zoo::coredata(time_serie) <- res
    
  } 
  
  return(time_serie)
  
}