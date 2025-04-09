qc_after_hmg <- function(target_data, mm_apply_qc = 0, apply_qc) {
  
  cubic_diff <- NULL
  new_cubic_diff <- NULL
  oro_val <- NULL
  
  # If lmt is NULL or NA, return data unchanged
  if (is.null(apply_qc) || is.na(apply_qc)) return(target_data)
  
  # Apply cubic root transformation
  oro_hmg <- cbind(target_data$raw_time_series, target_data$hmg_time_series)
  colnames(oro_hmg) <- c("oro", "hmg")
  time_step <- stats::time(oro_hmg[oro_hmg[, "hmg"] > mm_apply_qc])
  
  oro_ts <- target_data$raw_time_series[time_step]^(1 / 3)
  hmg_ts <- target_data$hmg_time_series[time_step]^(1 / 3)
  
  # Compute cubic differences
  cubic_dff <- hmg_ts - oro_ts
  n_cubic_dff <- cubic_dff[cubic_dff >= apply_qc | cubic_dff <= -apply_qc]
  
  # If no values exceed the threshold, return original data
  if (nrow(n_cubic_dff) == 0) return(target_data)
  
  # Prepare adjustments for exceeded values
  dates_time_n <-
    data.frame(
      hmg_val = as.numeric(target_data$hmg_time_series[time_step][stats::time(n_cubic_dff)]),
      oro_val = as.numeric(target_data$raw_time_series[time_step][stats::time(n_cubic_dff)]),
      cubic_diff = as.numeric(n_cubic_dff),
      time = stats::time(n_cubic_dff)
    )
  
  dates_time_n <- transform(
    dates_time_n,
    new_cubic_diff = ifelse(cubic_diff > 0, apply_qc, -apply_qc)
  )
  
  # Adjust cubic differences and compute new homogenized values
  dates_time_n <- transform(
    dates_time_n,
    new_hmg_val = (new_cubic_diff + oro_val^(1/3))^(3)
  )
  
  # Convert to xts and round
  new_hmg_val <- xts::xts(dates_time_n$new_hmg_val, dates_time_n$time)
  new_hmg_val <- round(new_hmg_val, 1)
  
  # Update target_data with new homogenized values
  new_target_data <- target_data
  new_target_data$hmg_time_series[stats::time(new_hmg_val)] <- new_hmg_val
  
  return(new_target_data)
}