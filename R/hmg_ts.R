hmg_ts <- function(prec, sts, neibs_max = 8, neibs_min = 3, thres = 1e+6, cor_neibs = 0.5, cleaning = FALSE, perc_break = 7, wet_day = 0, window_c = 15, apply_qc = 1, mm_apply_qc = 0, ncpu = 1){
  
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