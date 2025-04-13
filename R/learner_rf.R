#' Function of the random forest model (rf)

#' @description The function of the random forest model (rf) is used to predict wet/dry (classification) days and the amount of rain on wet days (regression). This function employs the hybrid model of Huerta et al. (2024).
#' @param ref matrix or data.frame containing the covariable data and precipitation value for each point location. This data is used to train (build) the model
#' @param can matrix or data.frame containing the variable data and precipitation value for one single point location. This data is used to execute the build model and get the predictions
#' @param covars vector of character elements containing the covariables names used in the model.
#' @details
#' The function should not be used directly, as its main purpose is to be passed as an argument in the core functions (qcPrec, gapFilling and gridPcp) of the package.
#' The core functions can handle whatever condition that avoids the execution of this function. Therefore, there is no need to test different situations on the data (empty values, among others). However, users who want to build their functions should test first before being passed as an argument in the core functions (you can check the other models (learner) functions available in the package).
#' This function encapsulates two types of models: classification and regression, and its output is a single numeric vector of three elements: probability of wet day, amount of wet day, and uncertainty of the amount of wet day.
#' To use this function, the randomForest package should have been installed previously
#' @references Huerta, A., Serrano-Notivoli, R., & Brönnimann, S. (2024). SC-PREC4SA: A serially complete daily precipitation dataset for South America. https://doi.org/10.31223/X57D8R
#' @export
#' @importFrom randomForest randomForest
#' @examples
#' \dontrun{
#' set.seed(123)
#' # creating random data (three predictors)
#' lon = rnorm(50,0,1)
#' lat = rnorm(50,40,1)
#' dcoast = rnorm(50,200,50)
#' prec = rnorm(1*50, mean = 1.2, sd = 6)
#' prec[prec<0] <- 0
#' 
#' # precipitation column should be call as "val"
#' data_full <- data.frame(lon = lon, lat = lat, dcoast = dcoast, val = prec)
#' 
#' # parameters for learner_rf
#' ref = data_full[-1, ]
#' can = data_full[1, ]
#' covars = c("lon", "lat", "dcoast")
#' 
#' learner_rf(ref = ref, can = can, covars = covars)
#' 
#' # case when "prec" is full wet
#' prec = rnorm(1*50, mean = 1.2, sd = 6)
#' prec[prec<0] <- 1
#' 
#' data_full <- data.frame(lon = lon, lat = lat, dcoast = dcoast, val = prec)
#' 
#' # parameters for learner_rf
#' ref = data_full[-1, ]
#' can = data_full[1, ]
#' covars = c("lon", "lat", "dcoast")
#' 
#' learner_rf(ref = ref, can = can, covars = covars)
#' }
#' 

learner_rf <- function(ref, can, covars) {
  
  #####################################  
  # probability of ocurrence prediction
  #####################################

  rr <- as.data.frame(ref)
  rr$val[rr$val > 0] <- 1
  rr$val <- factor(rr$val, levels = c(0, 1))
  
  # model
  f <- as.formula(
    paste("val ~ ", paste(covars, collapse = " + "), collapse = " ")
  )
  
  # prediction
  ## condition because original model can not deal with an unique factor
  if (length(unique(rr$val)) != 1) {
    
    set.seed(123)
    fmtb <- randomForest::randomForest(f,
                                       data = rr[, c("val", covars)],
                                       type = "classification")
    pb <- predict(fmtb, as.data.frame(can[, covars]), type = "prob")
    pb <- data.frame(pb)
    pb <- pb[, "X1"]
    
  } else {
    
    pb <- as.numeric(
      levels(unique(rr$val))
    )[unique(rr$val)]
    
  }
  
  #####################################  
  # amount prediction
  #####################################
  
  # rescaling
  rr <- as.data.frame(ref)
  MINc <- min(rr$val) -
    (as.numeric(quantile(rr$val, 0.50)) - as.numeric(quantile(rr$val, 0.25)))
  MINc <- ifelse(MINc < 0, 0, MINc)
  MAXc <- max(rr$val) +
    (as.numeric(quantile(rr$val, 0.75)) - as.numeric(quantile(rr$val, 0.50)))
  RANGE <- as.numeric(MAXc - MINc)
  rr$val <- (rr$val - MINc) / RANGE
  
  # model
  set.seed(123)
  fmt <- randomForest::randomForest(f,
                                    data = rr[, c("val", covars)],
                                    type = "regression")
  
  # prediction
  p <- predict(fmt, newdata = as.data.frame(can[, covars]), type = "response")
  p <- (p * RANGE) + MINc
  
  # error calculation
  e <- sqrt(
    sum(
      (rr$val - predict(fmt,
                        newdata = rr[, c("val", covars)],
                        type = "response")) ^ 2
    ) /
      (length(rr$val) - length(covars))
  )
  e <- (e * RANGE) + MINc
  
  #####################################  
  # out
  #####################################
  
  out <- as.numeric(c(pb, p, e))
  return(out)
  
}
