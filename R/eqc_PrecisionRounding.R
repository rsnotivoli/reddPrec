#' Precision and rounding patterns of daily precipitation

#' @description The function determine the level (0, 1 or 2) of the precision and rounding patterns test
#' @param xts_obj xts of a single time series
#' @param lmn_yday numeric value of the minimum number of days to be considered a complete year. The default value is 365 * 80 / 100 days
#' @param ths list. List of parameters to define in which level would be the daily precipitation.
#' @export
#' @details
#' Precision and rounding patterns depict inconsistencies in the frequency of decimal values in the time series. As there is no absolute correct frequency of decimals, we decided to measure how similar the decimal patterns are in the time series. A decimal pattern is interpreted as the list of unique decimal values observed, sorted in descending order. In this way, the decimal pattern for each year is computed first, followed by the selection of the most dominating pattern (mode). Based on how much (in percentage) this dominating pattern represents the time series, it is defined:
#' Level 0: coherent precision and rounding pattern (similar decimal pattern in more (or equal) than 70% of the time series).
#' Level 1: a similar decimal pattern in less than 70% but more than (or equal) 50% of the time series.
#' Level 2: different decimal patterns (no dominant pattern).
#' Argument ths present the above default thresholds for the levels definition.
#' @references Hunziker, S., Gubler, S., Calle, J., Moreno, I., Andrade, M., Velarde, F., ... & Brönnimann, S. (2017). Identifying, attributing, and overcoming common data quality issues of manned station observations. https://doi.org/10.1002/joc.5037
#' @references Huerta, A., Serrano-Notivoli, R., & Brönnimann, S. (2024). SC-PREC4SA: A serially complete daily precipitation dataset for South America. https://doi.org/10.31223/X57D8R
#' @examples
#' \dontrun{
#' set.seed(123)
#' 
#' # creating fake daily precipitation data
#' dates_t <- seq(as.Date("1980-01-01"), as.Date("2015-12-31"), by = "day")
#' prec <- round(rnorm(length(dates_t), mean = 1.2, sd = 6), 1)
#' prec[prec<0] <- 0
#' xts_obj <- xts::xts(prec, dates_t)
#' names(xts_obj) <- "prec"
#' 
#' # Comparison with visual inspection
#' eqc_Plot(xts_obj)
#' eqc_PrecisionRounding(xts_obj)
#' 
#' # it also work if there is some empty data (but not if all is NA)
#' xts_obj["1990/2010"] <- NA
#' eqc_Plot(xts_obj)
#' eqc_PrecisionRounding(xts_obj)
#' 
#' # creating some rounding values
#' xts_obj["1980/1990"] <- round(xts_obj["1980/1990"], 0) 
#' eqc_Plot(xts_obj)
#' eqc_PrecisionRounding(xts_obj)
#' 
#' }
#' 


eqc_PrecisionRounding <- function(xts_obj,
                                  lmn_yday = 365 * 80 / 100,
                                  ths = list(lv0 = c("percent" = 70),
                                             lv1to2 = c("percent" = 50))) {
  
  if (is.null(xts_obj) || all(is.na(xts_obj))) {
    stop("Error: The xts object is NULL or completely NA.")
  }
  
  # Get the decimal patterns from the xts object
  out <- get_dec_patterns(xts_obj = xts_obj, lmn_yday = lmn_yday)
  
  # Summarize the pattern shapes
  pattern_shape_types <- summary(out$pattern_shape)
  
  # Get the mode (most frequent value) of the pattern shapes
  pattern_shape_mode <- max(pattern_shape_types)
  
  # Determine the response based on the mode frequency
  # Level 0
  if (pattern_shape_mode >= (nrow(out) * ths$lv0["percent"] / 100)) {
    response <- 0

  # Level 1
  } else if (pattern_shape_mode >= (nrow(out) * ths$lv1to2["percent"] / 100)) {
    response <- 1

  # Level 2
  } else {
    response <- 2
  }
  
  return(response)

}
