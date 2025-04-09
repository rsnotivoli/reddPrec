
get_ndec <- function(xts_obj, lmn_yday, toPlot = FALSE) {
  
  # Calculate the absolute values and fractional part
  w <- abs(xts_obj)
  w <- (round(w, 1) - as.integer(round(w, 1)))
  w <- as.integer(w * 10)
  
  # Create a data frame with 'year' and 'dec' (factor of decimals)
  out_df <- data.frame(
    year = as.numeric(format(stats::time(xts_obj), "%Y")),
    dec = factor(
      w,
      levels = seq(9, 0, -1),
      labels = paste("x.", seq(9, 0, -1), sep = "")
    )
  )
  
  if(isTRUE(toPlot)){
    
    return(out_df)
    
  } else {
    
    # Remove rows with missing values
    out_df <- out_df[stats::complete.cases(out_df), ]
    
    # Reshape the data to wide format, counting occurrences of each decimal level
    out <- reshape2::dcast(out_df, year ~ dec,
                           value.var = "dec",
                           fun.aggregate = length)
    
    # Calculate 'size' as the row sum of the decimal counts (excluding 'year' column)
    if (ncol(out) > 2) {
      out$size <- rowSums(out[,-1])
    } else {
      out$size <- out[, 2]
    }
    
    # Filter rows where 'size' is greater than the specified threshold
    out <- out[out$size > lmn_yday, ]
    
    return(out)
    
  }

}