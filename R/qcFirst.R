
#' This function iterates the Quality Control by single time steps,
#' it can be parallelized

#' @param x vector of precipitation (original) of a single time step and an extra figure with the previous numebr of detected supects
#' @param distanc matrix of distances with the names of the neighbouring stations
#' @param it integer number of iteration
#' @param sts data.frame with the stations info
#' @param neibs integer number of nearest neighbours to use

.qcFirst <- function(x, distanc, it, sts, neibs) {
  
  detect <- x[length(x)] #number of suspects
  x <- x[-c(length(x))] #data
  nams <- as.character(sts$ID) #names
  
  w <- which(!is.na(x))
  if(length(w) == 0) { #if no data, ends
    return(c(x, 0))
    } else{
    
      #subset dataset based on available data
      sub <- x[w]
      nams <- nams[w]
      dd <- distanc[match(nams, rownames(distanc)),]
      
   if (it > 1 & detect == 0) #if no more suspects ends
      return(c(x, 0))
   else {
      if (max(sub, na.rm = T) == 0) {
        return(c(x, 0))
      }
      else {
        
        #start evaluating data
        code <- rep(NA, length(sub))
        for (h in 1:length(sub)) {
          #candidate
          can <- sub[h]
          #set nearest observations
          ww <- match(nams[-c(h)],dd[h,])
          ww <- sort(ww)[1:(neibs+1)]
          dst <- dd[h,][ww]
          ref <- sub[match(dst, nams)]
          if (max(ref) == 0 & can > 0) {
                code[h] <- 1
          }
          else if (min(ref) > 0 & can == 0) {
                code[h] <- 2
            }
            else {
              if(max(ref) == 0){
                pb <- 0
                p <- 0
              }
              else{
                sts_can <- sts[which(nams[h] == sts$ID), ]
                alts <- sts$ALT[match(dst, sts$ID)]
                lats <- sts$LAT[match(dst, sts$ID)]
                lons <- sts$LON[match(dst, sts$ID)]
                
                #probability of ocurrence prediction
                b <- ref
                b[b > 0] = 1
                fmtb <- suppressWarnings(glm(b ~ alts +
                                               lats + 
                                               lons, 
                                               family = binomial()))
                newdata = data.frame(alts = sts$ALT[which(nams[h] == sts$ID)], 
                                     lats = sts$LAT[which(nams[h] == sts$ID)], 
                                     lons = sts$LON[which(nams[h] == sts$ID)])
                pb <- predict(fmtb, newdata = newdata, 
                              type = "response")
                if (pb < 0.001 & pb > 0) 
                  pb <- 0.001
                pb <- round(pb, 3)
                
                #amount prediction
                #rescaling
                mini <- min(ref)/2
                maxi <- max(ref) + (max(ref) - min(ref))
                yr <- as.numeric((ref - mini)/(maxi - mini))
                fmt <- suppressWarnings(glm(yr ~ alts + 
                                            lats + 
                                            lons, 
                                            family = quasibinomial()))
                p <- predict(fmt, newdata = newdata, 
                             type = "response")
                if (p < 0.001 & p > 0) 
                  p <- 0.001
                p <- round((p * (maxi - mini)) + mini, 3)
              }
              
                #evaluating outliers
                if (can == 0 & pb > 0.5) {
                  if ((max((can + 10)/(p + 10), 
                           (p + 10)/(can + 10))) > 10) {
                    code[h] <- 3
                  }
                }
                if (can > 0) {
                  if ((max((can + 10)/(p + 10), 
                           (p + 10)/(can + 10))) > 10) {
                    code[h] <- 3
                  }
                }
                
                #evaluating suspect dry
                if (can == 0 & pb > 0.99 & p > 50) {
                  code[h] <- 4
                }
                
                #evaluating suspect wet
                if (can > 50 & pb < 0.01 & p < 1) {
                  code[h] <- 5
                }
              }
            }

          #encoding detected
          ww = which(!is.na(code))
          if (length(ww) > 0) {
            x[match(unique(nams[ww]), colnames(prec))] <- NA
            detect <- length(ww)
          }
          else detect <- 0
          return(c(x, detect))
        }
    }
  }
}
