
#' This function applies a final Quality Control by single time steps

#' @param prec matrix or data.frame containing the original precipitation data. Each column represents one station. The names of columns have to be names of the stations.
#' @param sts matrix or data.frame. A column "ID" (unique ID of stations) is required. The rest of the columns (all of them) will act as predictors of the model.
#' @param crs character. Coordinates system in EPSG format (e.g.: "EPSG:4326").
#' @param coords vector of two character elements. Names of the fields in "sts" containing longitude and latitude.
#' @param coords_as_preds logical. If TRUE (default), "coords" are also taken as predictors.
#' @param neibs integer. Number of nearest neighbors to use.
#' @param thres numeric. Maximum radius (in km) where neighboring stations will be searched. NA value uses the whole spatial domain.
#' @param qc vector of strings with the QC criteria to apply. Default is "all". See details.
#' @param qc3 numeric. Indicates the threshold (number of times higher or lower) from which a observation, in comparison with its estimate, should be deleted. Default is 10.
#' @param qc4 numeric vector of length 2. Thresholds of wet probability (0 to 1) and magnitude (in the units of input precipitation data) from which a observation of value zero, in comparison with its estimate, should be deleted. Default is c(0.99, 5). 
#' @param qc5 numeric vector of length 2. Thresholds of dry probability (0 to 1) and magnitude (in the units of input precipitation data) from which a observation higher than a specific value (also in the original units), in comparison with its estimate, should be deleted. Default is c(0.01, 0.1, 5). 
#' @importFrom terra distance relate
#' @noRd
  

qcLast <- function(x, y, sts, neibs, coords, crs, coords_as_preds = TRUE, qc = 'all', qc3 = 10, qc4 = c(0.99, 5), qc5 = c(0.01, 0.1, 5), thres = NA) {
  
  pb <- p <- NA
  
  n <- names(sts)
  if(!coords_as_preds) n <- n[-match(coords,names(sts))]
  if(neibs<length(n)*2) 
    stop('The number of neighbors must be at least the double of predictors')
  
  n <- names(sts)
  if(!coords_as_preds) n <- n[-match(coords,names(sts))]
  covars <- paste(n, collapse='+') # predictors
  
  if(qc == "all") qc <- as.character(paste0(1:5))
  
  wy <- which(!is.na(y)) # available on original
  wx <- which(!is.na(x)) # available on cleaned
  
  cc <- y*NA # code of removed data
  
  # if no values on original dataset, 
    # or number of cleaned data is lower than neibs, ends.
  if(length(wy) == 0 | length(wx) < neibs) { 
    return(c(y,cc))
  } else {
    
    #subset dataset based on available data
    cle <- data.frame(sts[wx,],val = x[wx])
    cle <- vect(cle, geom = coords,crs = crs, keepgeom = TRUE)
    ori <- data.frame(sts[wy,],val = y[wy])
    ori <- vect(ori, geom = coords,crs = crs, keepgeom = TRUE)
    
    
    if (max(y, na.rm = T) == 0) {
      return(c(y,cc))
    }
      else {
        #start evaluating data
        code <- rep(NA, length(y))
        for (h in 1:length(ori)) {
          can <- ori[h,]
          # checks if candidate is within reference
          ws <- which(relate(can, cle, "intersects"))
          # if it is, pull it out
          if(length(ws)>0) ref <- cle[-ws] else ref <- cle
          # cleaned data act as neighboring observations
          dd <- distance(can,ref)/1000
          if(!is.na(thres)){ 
            dd <- dd[dd<thres]
          }
          if(length(dd)<neibs){
            # message(paste0("Not enough observations within radius"))
            pb <- p <- NA
          } else{
          ref <- ref[match(sort(dd)[1:neibs],dd)] # just the number of neibs
            
          #checks for codes
          if (max(ref$val) == 0 & can$val > 0) {
            if(grep("1",qc)>0) code[h] <- 1
          } else if (min(ref$val) > 0 & can$val == 0) {
            if(grep("2",qc)>0) code[h] <- 2
          } else {
            if(max(ref$val) == 0){
              pb <- 0
              p <- 0
            } else{
                # probability of ocurrence prediction
                rr <- as.data.frame(ref)
                rr$val[rr$val > 0] <- 1
                
                
                f <- as.formula(paste0('val ~ ',covars))
                fmtb <- suppressWarnings(
                  glm(f,family = binomial(),data = rr)
                )
                
                pb <- round(predict(fmtb, newdata = as.data.frame(can), 
                                    type = "response"),3)
                
                #amount prediction
                #rescaling
                rr <- as.data.frame(ref)
                MINc <- min(rr$val) - (as.numeric(quantile(rr$val, 0.50)) - as.numeric(quantile(rr$val, 0.25)))
                MINc <- ifelse(MINc<0,0,MINc)
                MAXc <- max(rr$val) + (as.numeric(quantile(rr$val, 0.75))-as.numeric(quantile(rr$val, 0.50)))
                RANGE <- as.numeric(MAXc - MINc)
                rr$val <- (rr$val - MINc) / RANGE
                
                fmt <- suppressWarnings(
                  glm(f,family = quasibinomial(),data = rr)
                )
                p <- predict(fmt, newdata = as.data.frame(can),type = "response")
                p <- round((p * RANGE) + MINc, 3)
              }
            }
          }
            if(is.na(pb) | is.na(p)){
              code[h] <- NA
            } else{
              #evaluating outliers
              if (can$val == 0 & pb > 0.5) {
                if(grep("3",qc)>0){
                  if ((max((can$val + 0.1)/(p + 0.1), 
                           (p + 0.1)/(can$val + 0.1))) > qc3) {
                    code[h] <- 3
                  }
                }
              }
              if (can$val > 0) {
                if(grep("3",qc)>0){
                  if ((max((can$val + 0.1)/(p + 0.1), # 0.1 avoids problems with zeros
                           (p + 0.1)/(can$val + 0.1))) > qc3) {
                    code[h] <- 3
                  }
                }
              }
              
              #evaluating suspect dry
              if (can$val == 0 & pb > qc4[1] & p > qc4[2]) {
                if(grep("4",qc)>0) code[h] <- 4
              }
              
              #evaluating suspect wet
              if (can$val > qc5[3] & pb < qc5[1] & p < qc5[2]) {
                if(grep("5",qc)>0) code[h] <- 5
              }
          }
        } # end of calculations for all observations
        
        #encoding detected
        ww <- which(!is.na(code))
        if (length(ww) > 0) {
          y[wy][ww] <- NA
          cc[wy][ww] <- code[ww]
        }
        return(c(y,cc))
        
      } # end of maximum is zero condition
    
  } # end of number of neibs condition
  
} # end of function
