#' Weekly cycle test data frame
#' 
#' @description The function get the results of the statistical test of the wet day fraction within the weekly cycle test
#' @param xts_obj xts of a single time series
#' @noRd
#' @references Hunziker, S., Gubler, S., Calle, J., Moreno, I., Andrade, M., Velarde, F., ... & Brönnimann, S. (2017). Identifying, attributing, and overcoming common data quality issues of manned station observations. https://doi.org/10.1002/joc.5037
#' @references Huerta, A., Serrano-Notivoli, R., & Brönnimann, S. (2024). SC-PREC4SA: A serially complete daily precipitation dataset for South America. https://doi.org/10.31223/X57D8R
#' 

get_wd_fraction <- function(xts_obj) {
  week <- NULL
  original_locale <- Sys.getlocale("LC_TIME")  # Save the current locale
  Sys.setlocale("LC_TIME", "C")  # Set locale to English
  
  out_df <- data.frame(
    value = as.numeric(xts_obj),
    week = base::weekdays(stats::time(xts_obj))
  )
  
  out_df_wd <- out_df[out_df$value >= 0.1, ]
  
  length_values <- stats::aggregate(
    value ~ week, data = out_df,
    function(x) length(x[!is.na(x)]), na.action = NULL
  )
  
  if (all(is.na(out_df_wd$value))) {
    lenght_wd <- length_values
  } else {
    lenght_wd <- stats::aggregate(
      value ~ week, data = out_df_wd,
      function(x) length(x[!is.na(x)]), na.action = NULL
    )
  }
  
  out_df <- merge(length_values, lenght_wd, by = "week")
  colnames(out_df) <- c("week", "count", "count_wd")
  
  out_df <- transform(
    out_df,
    week = factor(
      week,
      levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"),
      labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"),
      ordered = TRUE
    )
  )
  
  out_df$frac_wd <- out_df$count_wd / out_df$count
  
  if (is.nan(sum(out_df$frac_wd)) | is.na(sum(out_df$frac_wd))) {
    out_df$bin_test <- factor(NA, levels = c("Rejected Ho", "No Rejected Ho"))
  } else {
    for (i in 1:nrow(out_df)) {
      test_bt <- stats::binom.test(
        out_df$count_wd[i], out_df$count[i],
        p = sum(out_df$count_wd) / sum(out_df$count),
        alternative = "two.sided", conf.level = 0.95
      )
      
      out_df$bin_test[i] <- ifelse(
        !((test_bt$conf.int[1] < sum(out_df$count_wd) / sum(out_df$count)) &
            (sum(out_df$count_wd) / sum(out_df$count) < test_bt$conf.int[2])),
        "Rejected Ho", "No Rejected Ho"
      )
    }
    
    out_df$bin_test <- factor(
      out_df$bin_test, levels = c("Rejected Ho", "No Rejected Ho")
    )
  }
  
  out_df <- out_df[order(out_df$week), ]
  rownames(out_df) <- NULL
  out_df$frac_wd <- round(out_df$frac_wd * 100, 1)
  
  Sys.setlocale("LC_TIME", original_locale)  # Reset to the original locale
  return(out_df)
}