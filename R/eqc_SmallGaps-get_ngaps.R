#' Small gaps test data frame
#' 
#' @description The function get the results of the small gaps test
#' @param xts_obj xts of a single time series
#' @param lmn_yday numeric value of the minimum number of days to be considered a complete year. The default value is 365 * 80 / 100 days
#' @noRd
#' @references Hunziker, S., Gubler, S., Calle, J., Moreno, I., Andrade, M., Velarde, F., ... & Brönnimann, S. (2017). Identifying, attributing, and overcoming common data quality issues of manned station observations. https://doi.org/10.1002/joc.5037
#' @references Huerta, A., Serrano-Notivoli, R., & Brönnimann, S. (2024). SC-PREC4SA: A serially complete daily precipitation dataset for South America. https://doi.org/10.31223/X57D8R
#' 

get_ngaps <- function(xts_obj,
                      lmn_yday) {
  
  # Calculate the total number of non-missing values per year
  size_y <- xts::apply.yearly(xts_obj, function(x) sum(!is.na(x)))
  
  # Count values within specific ranges per year
  var_01 <- xts::apply.yearly(
    xts_obj,
    function(x) sum(!is.na(x[x < 1 & x > 0]))
  )
  var_12 <- xts::apply.yearly(
    xts_obj,
    function(x) sum(!is.na(x[x < 2 & x > 1]))
  )
  var_23 <- xts::apply.yearly(
    xts_obj,
    function(x) sum(!is.na(x[x < 3 & x > 2]))
  )
  var_34 <- xts::apply.yearly(
    xts_obj,
    function(x) sum(!is.na(x[x < 4 & x > 3]))
  )
  var_45 <- xts::apply.yearly(
    xts_obj,
    function(x) sum(!is.na(x[x < 5 & x > 4]))
  )
  
  # Combine all variables into a data frame
  var_df <- data.frame(
    year = as.numeric(format(stats::time(var_01), "%Y")),
    var01 = as.numeric(var_01),
    var12 = as.numeric(var_12),
    var23 = as.numeric(var_23),
    var34 = as.numeric(var_34),
    var45 = as.numeric(var_45),
    year_size = as.numeric(size_y)
  )
  
  # Filter rows based on the size of the year
  var_df <- var_df[var_df$year_size > lmn_yday, ]
  rownames(var_df) <- NULL
  
  return(var_df)
}