hmg_indices_ts <- function(target_data) {
  
  # Define month groups
  aprsep_months <- c("04", "05", "06", "07", "08", "09")
  octmar_months <- c("01", "02", "03", "10", "11", "12")
  
  # Helper function to apply yearly functions
  apply_yearly_func <- function(data, fun) {
    lapply(data, function(z) xts::apply.yearly(z, fun))
  }
  
  # Split data by months
  months <- format(stats::time(target_data), "%m")
  data_aprsep <- target_data[months %in% aprsep_months]
  data_octmar <- target_data[months %in% octmar_months]
 
  # Compute indices
  indices <- list(
    prcptot_annual = apply_yearly_func(target_data, function(x) sum(x, na.rm = TRUE)),
    prcptot_aprsep = apply_yearly_func(data_aprsep, function(x) sum(x, na.rm = TRUE)),
    prcptot_octmar = apply_yearly_func(data_octmar, function(x) sum(x, na.rm = TRUE)),
    r1mm_annual = apply_yearly_func(target_data, function(x) sum(x > 0.1, na.rm = TRUE)),
    r1mm_aprsep = apply_yearly_func(data_aprsep, function(x) sum(x > 0.1, na.rm = TRUE)),
    r1mm_octmar = apply_yearly_func(data_octmar, function(x) sum(x > 0.1, na.rm = TRUE))
  ) 
  
  # Convert lists to xts objects
  indices <- lapply(indices, function(lst) do.call(cbind, lst))
  
  return(list(raw = target_data,
              indices = indices))
  
}