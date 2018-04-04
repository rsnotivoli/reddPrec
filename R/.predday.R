

.predday <- function(x, dat, distanc, points, sts, datess, neibs, intermediate) {
  
  #check input data
  x[is.infinite(x)] <- NA
  
  nams <- as.character(sts$ID) #names
  at <- as.numeric(sts$ALT)#altitudes
  ln <- as.numeric(sts$LON)#longitudes
  lt <- as.numeric(sts$LAT)#latitudes
  
  output <- paste('./gridded/',as.character(datess[dat]),'.txt',sep='')
  #set output
  wd_pred <- rep(NA, nrow(points))
  pred <- rep(NA, nrow(points))
  err <- rep(NA, nrow(points))
  
  w <- which(!is.na(x))
  if(length(w) == 0) { #if no data, ends
    if(intermediate){
      write.table(data.frame(ID = points$ID,
                             pred = rep(NA, nrow(points)),
                             err = rep(NA, nrow(points))), output, quote =F, 
                  sep = '\t', row.names = F, na = '')}
    return(rep(NA, nrow(points)))
    print(paste('No data on day',as.character(datess[dat])))
  } else if(length(w) < neibs){ #if less data than neighbours, ends
    if(intermediate){
      write.table(data.frame(ID = points$ID,
                             pred = rep(NA, nrow(points)),
                             err = rep(NA, nrow(points))), output, quote =F, 
                  sep = '\t', row.names = F, na = '')}
    return(rep(NA, nrow(points)))
    print(paste('Not enough neighbours on day',as.character(datess[dat])))
    } else{
      if(max(x, na.rm = T) == 0) { #if max data is 0, ends
        if(intermediate){
          write.table(data.frame(ID = points$ID,
                                 pred = rep(0, nrow(points)),
                                 err = rep(0, nrow(points))), output, quote =F, 
                      sep = '\t', row.names = F, na = '')}
        return(rep(0, nrow(points)))
      } else {
        
        #start evaluating data
        for (h in 1:nrow(points)) {
          #set nearest observations
          dst <- distanc[h,]
          ref <- x[match(dst, nams)]
          alts <- at[match(dst, nams)]
          lats <- lt[match(dst, nams)]
          lons <- ln[match(dst, nams)]
          wna <- sum(is.na(ref))
          if(wna > 0){
            wa <- which(is.na(ref))
            ref <- ref[-c(wa)]
            alts <- alts[-c(wa)]
            lats <- lats[-c(wa)]
            lons <- lons[-c(wa)]
            dst <- dst[match(names(ref), dst)]}
          #if less than 5 values, ends
          if(length(which(!is.na(ref))) < 5){
            pred[h] <- NA
            err[h] <- NA
          } else{
            if (max(ref) == 0) {
              wd_pred[h] <- 0
              pred[h] <- 0
              err[h] <- 0
            } else{
              sts_can <- points[h, ]
              
              
              #probability of ocurrence prediction
              b <- as.numeric(ref)
              b[b > 0] = 1
              fmtb <- suppressWarnings(glm(b ~ alts +
                                             lats + 
                                             lons, 
                                           family = binomial(), na.action = na.omit))
              newdata = data.frame(alts = as.numeric(sts_can$ALT), 
                                   lats = as.numeric(sts_can$LAT), 
                                   lons = as.numeric(sts_can$LON))
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
              e <- sqrt(sum((yr - predict(fmt, type = 'response')) ^ 2)/(length(yr) - 3))
              e <- round((e * (maxi - mini)) + mini, 3)
              if(e < 0.001 & e > 0)
                e  <- 0.001
              
              #evaluating estimate
              wd_pred[h] <- pb
              pred[h] <- p
              err[h] <- e
              if(pb <= 0.5) pred[h] <- 0
            }
          }
        }
        if(intermediate){
          write.table(data.frame(ID = points$ID,
                                 pred = as.numeric(pred), 
                                 err = as.numeric(err)), 
                      output, quote =F, 
                      sep = '\t', row.names = F, na = '')}
        return(pred)
      }
  }
}
