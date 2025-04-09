
hmg_build_ts <- function(target_nearby,
                         neibs_min,
                         cor_neibs,
                         prec) {
  
  # Extract data for the given stations (target + nearby)
  out_xts <- prec[, target_nearby, drop = FALSE]
  
  # Compute yearly totals
  out_xts_yearly <- lapply(
    out_xts,
    function(idx) {
      xts::apply.yearly(idx, sum, na.rm = TRUE)
    }
  )
  out_xts_yearly <- do.call(cbind, out_xts_yearly)
  
  # Correlation based on yearly totals (originally was at daily values)
  cor_test <- stats::cor(out_xts_yearly)[1, -1] > cor_neibs
  cor_test <- cor_test[cor_test == 1]
  
  n_stations <- sum(cor_test, na.rm = TRUE)
  
  # Get correlated stations
  # Relative
  if (n_stations > neibs_min) {
    
    res <- prec[, c(target_nearby[1], names(cor_test))]
    
  } else {
  # Absolute
    res <- prec[, c(target_nearby[1])]
    
  }
  
  # Return the final dataset
  return(res)
  
}