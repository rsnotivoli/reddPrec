#' Homogenization of daily precipitation time series

#' @description The function applies the homogenization procedure (detection and adjustment) for daily precipitation based on Huerta et al. (2024).
#' @param prec xts matrix containing the raw (quality-controlled and gap-filled) precipitation data. Each column represents one station. The names of columns must coincide with the names of the stations.
#' @param sts data.frame. A column "ID" (unique ID of stations), "LON" (decimal degree), and "LAT" (decimal degree) are required.
#' @param neibs_max integer. Number of maximum nearest neighbors to use. NA or Null value uses the whole number of stations 
#' @param neibs_min integer. Number of minimum nearest neighbors to use. This value represent the option to use the absolute or relative approach on the homogenization.
#' @param thres numeric. Maximum radius (in m) where neighboring stations will be searched. NA or Null value uses the whole spatial domain.
#' @param cor_neibs numeric. Minimum value of temporal correlation to define nearest neighbors.
#' @param cleaning logical. Set to FALSE as default. TRUE if the time series should be cleaned (trends and autocorrelation removed) before the detection test, as in Lund et al. (2023).
#' @param perc_break numeric. Value that define the percentage of time series that are statistically significant in the detection test in order to define a break point (year).
#' @param wet_day numeric. Value that define if the adjustment should be only performed on wet day (> 0 mm). Negative value mean that the adjustment will be also applied on zeros.
#' @param window_c integer. Window size of the application of the quantile matching adjustment.
#' @param apply_qc numeric. Maximum threshold (cubic difference) in which the adjusted data be considered corrected. Set to 1 as default. If the difference of the root cubic between the adjusted and raw data is above 1, the adjusted value will be corrected to be not above (or below) 1. 
#' @param mm_apply_qc numeric. Precipitation threshold in which apply_qc would be applied. Set to 0 as default, meaning that only the adjustment will be on wet days (> 0 mm)
#' @param ncpu number of processor cores used to parallel computing.
#' @export
#' @importFrom xts xts apply.yearly
#' @importFrom zoo coredata rollapply
#' @importFrom geosphere distHaversine
#' @importFrom Kendall MannKendall
#' @importFrom pracma detrend
#' @importFrom car durbinWatsonTest
#' @importFrom BreakPoints pettit man.whi stu SNHT Buishand_R
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @details
#' The homogenization procedure uses an automatic algorithm for both detection and adjustment without metadata information. In addition, relative and absolute approaches are combined for situations in which relative homogenization can not be performed. The absolute test, which has a lower power of detection than the relative tests, is thus intended as a backup test for when a relative test is hardly possible.
#' To ensure high confidence in breakpoint detection, a combination of different statistical tests and intercomparison of their results was used. Five univariate breakpoint tests were applied: Student's, Mann-Whitney, Buishand-R, Pettit, and the Standard Normal Homogeneity Test. Depending on the availability of nearby stations for a target time series: relative and absolute. For the relative approach, the algorithm searches for up to neibs_max well-correlated (> cor_neibs) nearby stations within a three radius. Later, the five tests are applied to difference series (target minus nearby) created with three different temporal aggregations. Finally, the breakpoint is set to a certain year if it is found in at least perc_break (%) of the number of difference time series that are significant (p-value < 0.05), using a tolerance of +/- 1 year. The absolute approach is used if the algorithm detects fewer than neibs_min nearby stations or none at all.
#' In the adjustment, it was adapted to the quantile-matching technique outlined in Squintu et al. (2018). It should be mentioned that this algorithm was created for temperature data; therefore, we made some changes to be used for precipitation. Dry values (wet_day parameter) cannot be corrected, and wet values were transformed twice (square root and log) before the algorithm execution to force a normal distribution. Based on this consideration, the correction was applied in two ways: relative and absolute. For the relative approach, the adjustment factor was computed using the target and nearby time series of the detection stage. It is assumed that the data after the break is correct; thus, the correction is backward. For the absolute case, the adjustment factor was computed using the target time series. This can be seen as an application of quantile mapping as there are no nearby stations.
#' Adjustment of daily precipitation can influence the extreme tails. Therefore, adjusted values can be set to not exceed a limit (apply_qc - difference of the root cubic) with the raw data. The apply_qc can be applied to all precipitation values or above a specific threshold: mm_apply_qc. This procedure still keeps the extreme adjustment while preventing the creation of extremely excessive values.
#' The output of the function is a list that contains: the detection results and the adjusted time series for each station.
#' @references Huerta, A., Serrano-Notivoli, R., & Brönnimann, S. (2024). SC-PREC4SA: A serially complete daily precipitation dataset for South America. https://doi.org/10.31223/X57D8R
#' @references Lund, R. B., Beaulieu, C., Killick, R., Lu, Q., & Shi, X. (2023). Good practices and common pitfalls in climate time series changepoint techniques: A review. Journal of Climate, 36(23), 8041-8057.
#' @references Squintu, A. A., van der Schrier, G., Brugnara, Y., & Klein Tank, A. (2018). Homogenization of daily ECA&D temperature series. International journal of climatology, 39(3), 1243-1261.
#' 


hmg_Ts <- function(prec, sts, neibs_max = 8, neibs_min = 3, thres = 1e+6, cor_neibs = 0.5, cleaning = FALSE, perc_break = 7, wet_day = 0, window_c = 15, apply_qc = 1, mm_apply_qc = 0, ncpu = 2){
  
  if (!xts::is.xts(prec)) {
    stop("Input data is not an xts object.")
  }
  
  if (any(is.na(prec))) {
    stop("Error: The xts object contains missing values (NA). Please handle missing data before proceeding.")
  }
  
  if(is.null(colnames(prec))){
    message('Guessed dataset names in same order as Stations ID')
    colnames(prec) <- sts$ID
  }
  
  m <- match(colnames(prec), sts$ID)
  if(length(which(is.na(m))) > 0) 
    stop('Stations ID do not coincide with dataset names')
  
  # reorder stations
  sts <- sts[m,]
  
  message(paste0('[',Sys.time(),'] -', " Homogenization"))
  
  registerDoParallel(cores=ncpu)
  
  j <- NULL
  a <- foreach(
    j = 1:ncol(prec),
    .export=c("hmg_nearby_ts", "hmg_build_ts", "hmg_indices_ts", "hmg_detection_ts", "hmg_adjustment_ts", "hmg_cleaning_trend_ar1", "qc_after_hmg")) %dopar% {
    
      target <- names(prec[, j])
      step_01 <- hmg_nearby_ts(target = target, sts = sts, neibs_max = neibs_max, thres = thres)
      step_02 <- hmg_build_ts(target_nearby = step_01, neibs_min = neibs_min, cor_neibs = cor_neibs, prec = prec)
      step_03 <- hmg_indices_ts(target_data = step_02)
      step_04 <- hmg_detection_ts(target_data = step_03, cleaning = cleaning, neibs_min = neibs_min, perc_break = perc_break)
      step_05 <- hmg_adjustment_ts(target_data = step_04, neibs_min = neibs_min, window_c = window_c, wet_day = wet_day)
      step_06 <- qc_after_hmg(target_data = step_05, apply_qc = apply_qc, mm_apply_qc = mm_apply_qc)
      step_06

  }
  
  a_adj_res <- do.call(cbind, lapply(a, `[[`, "hmg_time_series"))
  a_dec_res <- do.call(rbind, lapply(a, `[[`, "det_results"))
  # all(a_dec_res$ID == sts$ID) == all(colnames(prec) == colnames(a_adj_res))
  
  message(paste0('[',Sys.time(),'] -', " END"))
  
  return(
    list(adj_data = a_adj_res,
         det_data = a_dec_res)
  )
  
}