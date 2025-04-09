
hmg_nearby_ts <- function(target, sts, neibs_max = NULL, thres = NULL) {

  # Extract target coordinates
  target_xy <- sts[sts$ID == target, c("LON", "LAT"), drop = FALSE]
  
  # Ensure neibs_max and thres have valid default values
  neibs_max <- ifelse(is.na(neibs_max) | is.null(neibs_max), nrow(sts) - 1, neibs_max)
  thres <- ifelse(is.na(thres) | is.null(thres), Inf, thres)
  
  # Compute distances
  sts$distance <- geosphere::distHaversine(target_xy, sts[, c("LON", "LAT")])
  
  # Filter by threshold and sort
  filtered_sts <- sts[sts$distance < thres, ]
  filtered_sts <- filtered_sts[order(filtered_sts$distance), ]
  
  # Select nearest neighbors (ensuring it doesn't exceed available stations)
  filtered_sts <- filtered_sts[1:(neibs_max + 1), ]
  filtered_sts <- filtered_sts[stats::complete.cases(filtered_sts), ]
  
  filtered_sts <- as.character(filtered_sts$ID)

  return(filtered_sts)
}