#' Small gaps level of daily precipitation

#' @description The function determine the level (0, 1 or 2) of the small gaps test
#' @param xts_obj xts of a single time series
#' @param lmn_yday numeric value of the minimum number of days to be considered a complete year. The default value is 365 * 80 / 100 days
#' @param ths list. List of parameters to define in which level would be the daily precipitation.
#' @export
#' @details
#' Small gaps can be seen as unreported precipitation events that result in a gap or a frequency reduction in values below a specific threshold. To define the small gaps, we calculated the total count of values in five precipitation ranges from 0-1, 1-2, 2-3, 3-4, and 4-5 mm (not including the values in the limits) for each year. Therefore, considering the percentage of years with zero counts:
#' Level 0: no small gaps (0%; years with at least one value in any of the precipitation ranges).
#' Level 1: small gaps in at least 20% of years.
#' Level 2: small gaps in more than 20% of years
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
#' prec[prec<3 & prec>2] <- 0
#' xts_obj <- xts::xts(prec, dates_t)
#' names(xts_obj) <- "prec"
#' 
#' # Comparison with visual inspection
#' eqc_Plot(xts_obj)
#' eqc_SmallGaps(xts_obj)
#' 
#' # it also work if there is some empty data (but not if all is NA)
#' xts_obj["1990/2010"] <- NA
#' eqc_Plot(xts_obj)
#' eqc_SmallGaps(xts_obj)
#' 
#' }
#' 

eqc_SmallGaps <- function(xts_obj,
                          lmn_yday = 365 * 80 / 100,
                          ths = list(lv0 = c("percent" = 0),
                                     lv1to2 = c("percent" = 20))) {
  
  if (is.null(xts_obj) || all(is.na(xts_obj))) {
    stop("Error: The xts object is NULL or completely NA.")
  }
  
  # Get the gaps data frame using the get_ngaps function
  df_ngaps <- get_ngaps(xts_obj = xts_obj, lmn_yday = lmn_yday)
  
  # Extract relevant columns for processing
  df_ngaps_levels <- df_ngaps[, -match(c("year", "year_size"), colnames(df_ngaps))]
  
  # Replace all non-zero values with NA, then assign 1 to the remaining values
  df_ngaps_levels[df_ngaps_levels > 0] <- NA
  df_ngaps_levels[df_ngaps_levels >= 0] <- 1
  
  # Calculate the total number of gaps and their percentage
  ngaps <- max(colSums(df_ngaps_levels, na.rm = TRUE))
  percent_ngaps <- ngaps * 100 / nrow(df_ngaps_levels)
  
  # Determine response based on the percentage of gaps (levels)
  # Level 0
  if (percent_ngaps <= ths$lv0["percent"]) {
    response <- 0
  
  # Level 1
  } else if (percent_ngaps <= ths$lv1to2["percent"]) {
    response <- 1
  
  # Level 2
  } else {
    response <- 2
  }
  
  # Return the response level
  return(response)
}