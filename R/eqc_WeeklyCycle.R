#' Weekly cycle level of daily precipitation

#' @description The function determine the level (0, 1 or 2) of the weekly cycle test
#' @param xts_obj xts of a single time series
#' @param ths list. List of parameters to define in which level would be the daily precipitation.
#' @export
#' @details
#' Weekly cycles are characterized by the occurrence of wet days that significantly differ between the days of the week. To compute the weekly cycles, first, for each day of the week, the probability of precipitation is calculated by dividing the total number of wet days by the total counts of values. Later, the number of wet days is tested by a two-sided binomial test (95 % confidence level). Based on how many days were significant, it was defined:
#' Level 0: no atypical weekly cycle (similar probability between the days of the week).
#' Level 1: at least two days present an atypical probability (significant test).
#' Level 2: more than two days present an atypical probability (significant test) or one day presents an extremely different probability (more than 10%).
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
#' eqc_WeeklyCycle(xts_obj)
#' 
#' # it also work if there is some empty data (but not if all is NA)
#' xts_obj["1990/2010"] <- NA
#' eqc_Plot(xts_obj)
#' eqc_WeeklyCycle(xts_obj)
#' 
#' }
#' 

eqc_WeeklyCycle <- function(xts_obj,
                            ths = list(lv0 = c("n_days" = 0),
                                       lv1to2 = c("n_days" = 2, "percent" = 10))) {
  
  if (is.null(xts_obj) || all(is.na(xts_obj))) {
    stop("Error: The xts object is NULL or completely NA.")
  }
  
  wc_res <- get_wd_fraction(xts_obj = xts_obj)
  wc_res_rej <- wc_res[wc_res$bin_test == "Rejected Ho", ]
  wc_res_norej <- wc_res[wc_res$bin_test == "No Rejected Ho", ]
  
  # Special case with station X9446 in ES (Aragon)
  if (length(wc_res_norej$frac_wd) == 0) {
    
    diff_rejnorej <- Inf
    
  } else {
    
    diff_rejnorej <- sapply(
      wc_res_rej$frac_wd,
      function(idd) abs(idd - wc_res_norej$frac_wd)
    )
    
    diff_rejnorej <- ifelse(length(diff_rejnorej) < 1, 0, max(diff_rejnorej))
    
  }
  
  # Level 0
  if (nrow(wc_res_rej) <= ths$lv0["n_days"]) {
    response <- 0
  
  # Level 1
  } else if (nrow(wc_res_rej) <= ths$lv1to2["n_days"] & diff_rejnorej < ths$lv1to2["percent"]) {
    response <- 1
  
  # Level 2
  } else {
    response <- 2
  }
  
  return(response)

}
