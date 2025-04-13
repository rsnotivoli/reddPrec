#' Truncation smoothed line
#' 
#' @description The function get the smoothed line of the daily precipitation
#' @param xts_obj xts of a single time series
#' @noRd
#' @references Hunziker, S., Gubler, S., Calle, J., Moreno, I., Andrade, M., Velarde, F., ... & Brönnimann, S. (2017). Identifying, attributing, and overcoming common data quality issues of manned station observations. https://doi.org/10.1002/joc.5037
#' @references Huerta, A., Serrano-Notivoli, R., & Brönnimann, S. (2024). SC-PREC4SA: A serially complete daily precipitation dataset for South America. https://doi.org/10.31223/X57D8R
#'

get_ts_borderline <- function(xts_obj) {
  
  # Remove rows with missing values
  xts_obj1 <- xts_obj[stats::complete.cases(xts_obj)]
  
  # Apply rolling window to find the maximum value with a window size of 180
  border_filter <- zoo::rollapply(
    zoo::as.zoo(xts_obj1),
    width = 180,
    FUN = function(x) max(x, na.rm = TRUE),
    align = "center",
    partial = TRUE
  )
  
  # Apply rolling window again to find the maximum with a smaller window size of 90
  border_filter <- zoo::rollapply(
    border_filter,
    width = 90,
    FUN = function(x) max(x, na.rm = TRUE),
    align = "center",
    partial = TRUE
  )
  
  # Apply rolling window one more time to find the median with a window size of 180
  border_filter <- zoo::rollapply(
    border_filter,
    width = 180 * 1,
    FUN = function(x) stats::median(x, na.rm = TRUE),
    align = "center",
    partial = TRUE
  )
  
  # Round the result to the nearest 10 (this was done to avoid truncated values that differ in decimals)
  border_filter <- round(border_filter, digits = -1)
  
  return(border_filter)
  
}