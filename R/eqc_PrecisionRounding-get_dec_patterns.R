

get_dec_patterns <- function(xts_obj, lmn_yday) {
  
  # Get the decimal data using get_ndec function
  out <- get_ndec(xts_obj = xts_obj, lmn_yday = lmn_yday)
  
  # Prepare the data: exclude 'year' and 'size' columns
  out_dec <- out[, -match(c("year", "size"), colnames(out)), drop = FALSE]
  
  # Set values greater than 1 to 1, and values equal to 0 to NA
  out_dec[out_dec > 1] <- 1
  out_dec[out_dec == 0] <- NA
  
  # Create the 'pattern_shape' column by applying a function row-wise
  out$pattern_shape <- sapply(1:nrow(out_dec), function(idd) {
    
    # Multiply each row by its corresponding column number (ignoring 'x.' prefix)
    idd_d <- out_dec[idd, ] * as.numeric(gsub("x.", "", colnames(out_dec)))
    
    # Flatten and remove NA values
    idd_d <- unlist(idd_d)
    idd_d <- idd_d[!is.na(idd_d)]
    
    # Concatenate the remaining values into a pattern string
    paste0(idd_d, collapse = ".")
  })
  
  # Create the 'pattern_lenght' column by calculating the length of each pattern
  out$pattern_lenght <- sapply(1:nrow(out_dec), function(idd) {
    
    # Multiply each row by its corresponding column number (ignoring 'x.' prefix)
    idd_d <- out_dec[idd, ] * as.numeric(gsub("x.", "", colnames(out_dec)))
    
    # Flatten and remove NA values
    idd_d <- unlist(idd_d)
    idd_d <- idd_d[!is.na(idd_d)]
    
    # Return the length of the cleaned pattern
    length(idd_d)
  })
  
  # Convert 'pattern_shape' to a factor
  out$pattern_shape <- factor(out$pattern_shape)
  
  # Convert 'pattern_lenght' to a factor (not used, better only with pattern_shape)
  # out$pattern_lenght <- factor(out$pattern_lenght)
  
  return(out)
  
}