#' Truncation test data frame
#' 
#' @description The function get the results of the truncation test
#' @param xts_obj xts of a single time series
#' @importFrom reshape2 dcast
#' @noRd
#' @references Hunziker, S., Gubler, S., Calle, J., Moreno, I., Andrade, M., Velarde, F., ... & Brönnimann, S. (2017). Identifying, attributing, and overcoming common data quality issues of manned station observations. https://doi.org/10.1002/joc.5037
#' @references Huerta, A., Serrano-Notivoli, R., & Brönnimann, S. (2024). SC-PREC4SA: A serially complete daily precipitation dataset for South America. https://doi.org/10.31223/X57D8R
#'

get_lenght_borderline <- function(xts_obj) {
  
  # Get the borderline using the get_ts_borderline function
  border_line <- get_ts_borderline(xts_obj = xts_obj)
  
  # Remove missing values
  border_line_nona <- border_line[stats::complete.cases(border_line)]
  
  # Create a data frame with the borderline values and corresponding years
  border_line_df <- data.frame(
    value = factor(as.numeric(border_line_nona)),
    year = as.numeric(format(stats::time(border_line_nona), "%Y"))
  )
  
  # Reshape the data frame to get the count of each value by year
  border_line_df <- reshape2::dcast(
    border_line_df,
    year ~ value,
    fun.aggregate = length
  )
  
  # Add the size of each year (number of non-NA values)
  border_line_df$size <- as.numeric(
    xts::apply.yearly(
      border_line_nona,
      function(idd) sum(!is.na(idd))
    )
  )
  
  return(border_line_df)
  
}