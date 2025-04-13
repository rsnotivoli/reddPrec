#' Automatic enhanced quality control for daily precipitation time series 

#' @description The function determine the level (0, 1 or 2) of the enhanced quality control tests
#' @param prec xts matrix of precipitation time series
#' @param sts data.frame with metadata of the stations. A column "ID" (unique ID of stations) is required.
#' @param lmn_yday numeric. value of the minimum number of days to be considered a complete year. The default value is 365 * 80 / 100 days
#' @param ths_trc list. List of parameters to define in which level would be the daily precipitation for truncation.
#' @param ths_sgs list. List of parameters to define in which level would be the daily precipitation for small gaps.
#' @param ths_wcc list. List of parameters to define in which level would be the daily precipitation for weekly cycle.
#' @param ths_prp list. List of parameters to define in which level would be the daily precipitation for precision and rounding patterns.
#' @param ncpu integer. number of CPU threads to use for parallel computing.
#' @export
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @details
#' Truncation is when heavy precipitation episodes are truncated or noticeably reduced in frequency above a given threshold. Because there is no preceding algorithm for truncation, it is defined here as when the maximum boundary of a time series lasts for a set length of time (years). The maximum boundary is computed as the daily precipitation's maximum moving window value. Thus, based on the length of years:
#' Level 0: no truncation (a constant maximum boundary lasts at least 2 years).
#' Level 1: a constant maximum boundary lasts longer equal (or above) 3 years but less than 5 years
#' Level 2: a constant maximum boundary lasts more than 5 years.
#' Small gaps can be seen as unreported precipitation events that result in a gap or a frequency reduction in values below a specific threshold. To define the small gaps, we calculated the total count of values in five precipitation ranges from 0-1, 1-2, 2-3, 3-4, and 4-5 mm (not including the values in the limits) for each year. Therefore, considering the percentage of years with zero counts:
#' Level 0: no small gaps (0%; years with at least one value in any of the precipitation ranges).
#' Level 1: small gaps in at least 20% of years.
#' Level 2: small gaps in more than 20% of years
#' Weekly cycles are characterized by the occurrence of wet days that significantly differ between the days of the week. To compute the weekly cycles, first, for each day of the week, the probability of precipitation is calculated by dividing the total number of wet days by the total counts of values. Later, the number of wet days is tested by a two-sided binomial test (95 % confidence level). Based on how many days were significant, it was defined:
#' Level 0: no atypical weekly cycle (similar probability between the days of the week).
#' Level 1: at least two days present an atypical probability (significant test).
#' Level 2: more than two days present an atypical probability (significant test) or one day presents an extremely different probability (more than 10%).
#' Precision and rounding patterns depict inconsistencies in the frequency of decimal values in the time series. As there is no absolute correct frequency of decimals, we decided to measure how similar the decimal patterns are in the time series. A decimal pattern is interpreted as the list of unique decimal values observed, sorted in descending order. In this way, the decimal pattern for each year is computed first, followed by the selection of the most dominating pattern (mode). Based on how much (in percentage) this dominating pattern represents the time series, it is defined:
#' Level 0: coherent precision and rounding pattern (similar decimal pattern in more (or equal) than 70% of the time series).
#' Level 1: a similar decimal pattern in less than 70% but more than (or equal) 50% of the time series.
#' Level 2: different decimal patterns (no dominant pattern).
#' @references Hunziker, S., Gubler, S., Calle, J., Moreno, I., Andrade, M., Velarde, F., ... & Brönnimann, S. (2017). Identifying, attributing, and overcoming common data quality issues of manned station observations. https://doi.org/10.1002/joc.5037
#' @references Huerta, A., Serrano-Notivoli, R., & Brönnimann, S. (2024). SC-PREC4SA: A serially complete daily precipitation dataset for South America. https://doi.org/10.31223/X57D8R

eqc_Ts <- function(prec,
                   sts,
                   lmn_yday = 365 * 80 / 100,
                   ths_trc = list(lv0 = c("n_years" = 2), lv1to2 = c("n_years_l" = 3, "n_years_h" = 5)),
                   ths_sgs = list(lv0 = c("percent" = 0), lv1to2 = c("percent" = 20)), 
                   ths_wcc = list(lv0 = c("n_days" = 0), lv1to2 = c("n_days" = 2, "percent" = 10)),
                   ths_prp = list(lv0 = c("percent" = 70), lv1to2 = c("percent" = 50)),
                   ncpu = 1){
  
  if (is.null(prec) || all(is.na(prec))) {
    stop("Error: The xts object is NULL or completely NA.")
  }
  
  if(is.null(colnames(prec))){
    message('Guessed dataset names in same order as Stations ID')
    colnames(prec) <- sts$ID
  }
  
  message(paste0('[',Sys.time(),'] -', " Enhanced QC"))

  registerDoParallel(cores=ncpu)
  
  j <- NULL
  a <- foreach(
    j = 1:ncol(prec),
    .export=c("eqc_Truncation", "eqc_SmallGaps", "eqc_WeeklyCycle", "eqc_PrecisionRounding")) %dopar% {
      
      target <- prec[, j]
      
      eqc_trc <- eqc_Truncation(xts_obj = target, ths = ths_trc)
      eqc_sgs <- eqc_SmallGaps(xts_obj = target, ths = ths_sgs, lmn_yday = lmn_yday)
      eqc_wcc <- eqc_WeeklyCycle(xts_obj = target, ths = ths_wcc)
      eqc_prp <- eqc_PrecisionRounding(xts_obj = target, ths = ths_prp, lmn_yday = lmn_yday)
      
      data.frame(ID = names(target), Truncation = eqc_trc, Small_Gaps = eqc_sgs, Weekly_Cycle = eqc_wcc, Precision_Rounding = eqc_prp)
      
      }
  
  a <- do.call(rbind, a)
  
  message(paste0('[',Sys.time(),'] -', " END"))
  
  return(
    a  
  )
  
}