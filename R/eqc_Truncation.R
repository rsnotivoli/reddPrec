#' Truncation level of daily precipitation

#' @description The function determine the level (0, 1 or 2) of the truncation test
#' @param xts_obj xts of a single time series
#' @param ths list. List of parameters to define in which level would be the daily precipitation.
#' @export
#' @details
#' Truncation is when heavy precipitation episodes are truncated or noticeably reduced in frequency above a given threshold. Because there is no preceding algorithm for truncation, it is defined here as when the maximum boundary of a time series lasts for a set length of time (years). The maximum boundary is computed as the daily precipitation's maximum moving window value. Thus, based on the length of years:
#' Level 0: no truncation (a constant maximum boundary lasts at least 2 years).
#' Level 1: a constant maximum boundary lasts longer equal (or above) 3 years but less than 5 years
#' Level 2: a constant maximum boundary lasts more than 5 years.
#' Argument ths present the above default thresholds for the levels definition.
#' @references Hunziker, S., Gubler, S., Calle, J., Moreno, I., Andrade, M., Velarde, F., ... & Brönnimann, S. (2017). Identifying, attributing, and overcoming common data quality issues of manned station observations. https://doi.org/10.1002/joc.5037
#' @references Huerta, A., Serrano-Notivoli, R., & Brönnimann, S. (2024). SC-PREC4SA: A serially complete daily precipitation dataset for South America. https://doi.org/10.31223/X57D8R
#' @examples
#' \dontrun{
#' set.seed(123)
#' 
#' # creating fake daily precipitation data
#' dates_t <- seq(as.Date("1980-01-01"), as.Date("2015-12-31"), by = "day")
#' prec <- round(rnorm(length(dates_t), mean = 1.2, sd = 50), 1)
#' prec[prec<0] <- 0
#' prec[prec<3 & prec>2] <- 0
#' xts_obj <- xts::xts(prec, dates_t)
#' names(xts_obj) <- "prec"
#' 
#' # Comparison with visual inspection
#' eqc_Plot(xts_obj)
#' eqc_Truncation(xts_obj)
#' 
#' # it also work if there is some empty data (but not if all is NA)
#' xts_obj["1990/2010"] <- NA
#' eqc_Plot(xts_obj)
#' eqc_Truncation(xts_obj)
#' 
#' }
#' 

eqc_Truncation <- function(xts_obj,
                           ths = list(lv0 = c("n_years" = 2),
                                      lv1to2 = c("n_years_l" = 3, "n_years_h" = 5))) {
  
  if (is.null(xts_obj) || all(is.na(xts_obj))) {
    stop("Error: The xts object is NULL or completely NA.")
  }
  
  # Get the border line data using the get_lenght_borderline function
  df_border <- get_lenght_borderline(xts_obj = xts_obj)
  
  # Subtract the 'size' column from all other columns to flag differences
  df_flagg <- df_border[, -match(c("year", "size"), colnames(df_border))] -
    df_border$size
  
  # Set flagged values to NA
  df_flagg[df_flagg != 0] <- NA
  
  # Calculate the maximum length of contiguous runs (truncated lengths)
  trunc_lenght <- sapply(df_flagg, function(idd) max(rle(idd)$lengths))
  trunc_lenght <- max(trunc_lenght)
  
  # Determine the response based on the truncation length (years)
  # Level 0
  if (trunc_lenght <= ths$lv0["n_years"]) {
    response <- 0
  
  # Level 1
  } else if (trunc_lenght >= ths$lv1to2["n_years_l"] & trunc_lenght < ths$lv1to2["n_years_h"]) {
    response <- 1

  # Level 2
  } else {
    response <- 2
  }
  
  return(response)

}