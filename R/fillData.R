#' Estimates new precipitation values based on nearest neighbors.
#' 
#' @description This function calculates precipitation values of a single day.
#' @param x vector of observations in a single day.
#' @param sts matrix or data.frame containing the coordinates and predictors..
#' @param crs character. Coordinates system in EPSG format (e.g.: "EPSG:4326").
#' @param coords vector of two character elements. Names of the fields in "sts" containing longitude and latitude.
#' @param coords_as_preds logical. If TRUE (default), "coords" are also taken as predictors.
#' @param neibs integer. Number of nearest neighbors to use.
#' @param thres numeric. Maximum radius (in km) where neighboring stations will be searched. NA value uses the whole spatial domain.
#' @importFrom terra vect distance
#' @importFrom stats as.formula glm binomial predict quantile quasibinomial
#' @noRd

  fillData <-function(x, sts, neibs, coords, crs, coords_as_preds = TRUE, thres){
    
    n <- names(sts)
    if(!coords_as_preds) n <- n[-match(coords,names(sts))]
    if(neibs<length(n)*2) 
      stop('The number of neighbors must be at least the double of predictors')
    covars <- paste(n, collapse='+') # predictors
    
        #set output
    res <- matrix(NA, ncol = 6, nrow = length(x))
    colnames(res) <- c("obs","wd_pred", "raw_pred", 
                       "mod_pred", "err", "neibs")
    res[, 1] <- x
    
    w <- which(!is.na(x))
    if(length(w) == 0) { # if no data, ends
      return(res)
    } else if(length(w) < (neibs+1)){ # if less data than neighbours, ends
      return(res)
      } else{
        if(max(x, na.rm = T) == 0) { # if max value is 0, ends
          res[,2:5] <- 0
          return(res)
        } else {
          
          #subset dataset based on available data
          sub <- data.frame(sts,val = x)
          sub <- vect(sub, geom = coords, crs = crs, keepgeom = TRUE)
          
          ###################
          # evaluating data
          ###################
          for (h in 1:length(x)) {
            #candidate
            can <- sub[h]
            ref <- sub[-h]
            ref <- ref[!is.na(ref$val),]
            #set nearest observations
            dd <- terra::distance(can,ref)/1000
            if(!is.na(thres)){ 
              dd <- dd[dd<thres]
            }
            if(length(dd)<neibs){
              # message(paste0("Not enough observations within radius"))
              res[h,6] <- length(dd)
            } else{
              ref <- ref[match(sort(dd)[1:neibs],dd)]
              res[h,6] <- neibs
              if (max(ref$val) == 0) {
                res[h,2:5] <- 0
              } else if (sum(diff(ref$val))==0){
                res[h,2:5] <- c(1,ref$val[1],ref$val[1],0)
              } else{
                
                # probability of ocurrence prediction
                rr <- as.data.frame(ref)
                rr$val[rr$val > 0] <- 1
                
                
                f <- as.formula(paste0('val ~ ',covars))
                fmtb <- suppressWarnings(
                  glm(f,family = binomial(),data = rr)
                )
                
                pb <- round(predict(fmtb, newdata = as.data.frame(can), 
                                    type = "response"),2)
                
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
                p <- round((p * RANGE) + MINc, 2)
                
                # error calculation 
                e <- sqrt(sum((rr$val - predict(fmt, type = 'response')) ^ 2)/(length(rr$val) - length(n)))
                e <- round((e * RANGE) + MINc, 2)
                
                #evaluating estimate
                res[h, 2] <- pb
                res[h, 3] <- p
                if(pb <= 0.5) res[h, 4] <- 0 else res[h, 4] <- p
                res[h, 5] <- e
                
              } # end of is not zero
            } # end of enough neibs
          } # end of day
          
          return(res)
        } # end of max is zero
      } # end of not enough neibs
  } # EOF
  
  