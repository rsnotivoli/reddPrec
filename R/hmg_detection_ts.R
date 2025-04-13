hmg_detection_ts <- function(target_data, cleaning, neibs_min, perc_break){
  
  out <- lapply(names(target_data$indices), function(name) {
    
    # Extract variable and season from name
    name_split <- unlist(strsplit(name, "_"))
    variable <- name_split[1]
    season <- name_split[2]
    
    # Apply break detection
    data.frame(
      apply_break_detection(time_series = target_data$indices[[name]], cleaning = cleaning, neibs_min = neibs_min),
      variable = variable,
      season = season
    )
  })
  
  # Combine results and remove incomplete cases
  out <- do.call(rbind, out)
  out <- out[stats::complete.cases(out), ]
  
  # Extract break year
  break_year <- get_break_year(detection_results_output = out, perc_break = perc_break)
  
  return(list(
    original = target_data$raw,
    detection_test = out,
    break_year = break_year
  ))
  
}


apply_break_detection <- function(time_series, cleaning, neibs_min) {
  
  # Determine the number of stations (columns excluding the first one [target])
  n_stations <- ncol(time_series) - 1
  
  if(isTRUE(cleaning)) {
    
    # Relative
    if (n_stations > neibs_min) {
      
      # Multi-station processing: compute differences from reference (first column [target])
      time_serie_samples <- time_series
      time_serie_samples <- lapply(
        time_serie_samples[, -1],
        function(x) (time_serie_samples[, 1] - x)
      )
      time_serie_samples <- do.call(cbind, time_serie_samples)
      # Cleaning
      time_serie_samples <- lapply(
        time_serie_samples,
        function(x) hmg_cleaning_trend_ar1(time_serie = x)
      )
      time_serie_samples <- do.call(cbind, time_serie_samples)
      # Apply break detection
      time_serie_samples <- lapply(
        time_serie_samples,
        function(x) break_detection_tests(time_serie = x)
      )
      
      res <- do.call(rbind, time_serie_samples)

    # Absolute
    } else {
      
      # Single station processing
      time_serie_samples <- time_series[, 1]
      # Cleaning
      time_serie_samples <- hmg_cleaning_trend_ar1(time_serie = time_serie_samples)
      # Apply break detection
      res <- break_detection_tests(time_serie = time_serie_samples)
      
    }
     
  } else {
    
    # Relative
    if (n_stations > neibs_min) {
      
      # Multi-station processing: compute differences from reference (first column [target])
      time_serie_samples <- time_series
      time_serie_samples <- lapply(
        time_serie_samples[, -1],
        function(x) (time_serie_samples[, 1] - x)
      )
      time_serie_samples <- do.call(cbind, time_serie_samples)
      # Apply break detection
      time_serie_samples <- lapply(
        time_serie_samples,
        function(x) break_detection_tests(time_serie = x)
      )
      
      res <- do.call(rbind, time_serie_samples)
      
    # Absolute
    } else {
      
      # Single station processing
      time_serie_samples <- time_series[, 1]
      # Apply break detection
      res <- break_detection_tests(time_serie = time_serie_samples)

    }
    
  }
  
  return(res)
}

break_detection_tests <- function(time_serie, p_value = 0.05) {
  
  # Check if all values are the same
  if (identical(stats::var(as.numeric(time_serie)), 0)) {
    return(data.frame(
      test = "No test",
      breaks = NA,
      year_break = NA,
      p.value = NA,
      sig = NA
    ))
  }
  
  # List of break detection tests
  tests <- list(
    Pettitt = BreakPoints::pettit,
    Mann_Whitney_Wilcoxon = BreakPoints::man.whi,
    Studen_t = BreakPoints::stu,
    SNHT = function(serie) BreakPoints::SNHT(serie, simulations = 100),
    Buishand = function(serie) BreakPoints::Buishand_R(serie, simulations = 100)
  )
  
  # Apply each test and create a data.frame
  res <- do.call(rbind, lapply(names(tests), function(test_name) {
    
    test_result <- suppressWarnings(tests[[test_name]](serie = time_serie))
    data.frame(test_result, test = paste(test_name, "test"))
    
  }))
  
  # Extract break year
  res$year_break <- as.numeric(format(stats::time(time_serie)[res$breaks], "%Y"))
  
  # Significance column
  res$sig <- as.integer(res$p.value >= p_value)  # 1 if p >= p_value, 0 otherwise
  
  return(res[, c("test", "breaks", "year_break", "p.value", "sig")])
  
}

get_break_year <- function(detection_results_output, perc_break) {
  
  # Filter only significant breaks (sig == 0)
  res <- detection_results_output[detection_results_output$sig == 0, ]
  limit_n_test <- nrow(detection_results_output) * perc_break / 100
  
  if (nrow(res) == 0) return(NA)
  
  # Compute frequency of break occurrences within ±1 year range
  res$freq <- sapply(res$breaks, function(x) {
    length(res$breaks[res$breaks <= (x + 1) & res$breaks >= (x - 1)])
  })
  
  # Filter based on frequency threshold
  res <- res[res$freq > limit_n_test, ]
  # it is possible that after this point the data.frame does not have any row
  
  if (nrow(res) == 0) return(NA)
  
  # Determine the most frequent break year
  return(
    getmode(res$year_break[res$freq == max(res$freq, na.rm = TRUE)])
  )
}

getmode <- function(v) {
  
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
  
}