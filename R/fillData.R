

  fillData <-function(x, dat, datess, neibs, distanc = distanc, sts = sts, intermediate){
    
    nams <- as.character(sts$ID) #names
    dir.create('./days/',showWarnings = F)
    output <- paste('./days/',as.character(datess[dat]),'.txt',sep='')
    #set output
    wd_pred <- rep(NA, length(x))
    raw_pred <- rep(NA, length(x))
    mod_pred <- rep(NA, length(x))
    err <- rep(NA, length(x))
    
    w <- which(!is.na(x))
    if(length(w) == 0) { # if no data, ends
      if(intermediate){
      write.table(data.frame(ID = nams, 
                             obs = x, wd_pred = NA, 
                             raw_pred = NA, mod_pred = NA,
                             err = NA), output, quote =F, 
                  sep = '\t', row.names = F, na = '')}
      return(rep(NA, length(x)))
      message(paste('No data on day',as.character(datess[dat])))
    } else if(length(w) < neibs){ # if less data than neighbours, ends
      if(intermediate){
      write.table(data.frame(ID = nams, 
                             obs = x, wd_pred = NA, 
                             raw_pred = NA, mod_pred = NA,
                             err = NA), output, quote =F, 
                  sep = '\t', row.names = F, na = '')}
      return(rep(NA, length(x)))
      message(paste('Not enough neighbours on day',as.character(datess[dat])))
      } else{
        if(max(x, na.rm = T) == 0) { # if max value is 0, ends
          if(intermediate){
          write.table(data.frame(ID = nams, 
                                 obs = x, wd_pred = 0, 
                                 raw_pred = 0, mod_pred = 0,
                                 err = 0), output, quote =F, 
                      sep = '\t', row.names = F, na = '')}
          return(rep(0, length(x)))
        } else {
          
          ###################
          # evaluating data
          ###################
          for (h in 1:length(x)) {
            #set nearest observations
            dst <- distanc[h,]
            n <- nams[which(!is.na(x))]
            wn <- which(n == nams[h])
            if(length(wn) > 0) n <- n[-c(wn)]
            k <- match(dst, n)
            dst <- dst[-c(which(is.na(k)))] # neighbours with data
            if(length(dst) < neibs) next # ends if inssuficient neibs
            dst <- dst[1:neibs]
            ref <- x[match(dst, nams)]
            if (max(ref) == 0) {
              wd_pred[h] <- 0
              raw_pred[h] <- 0
              mod_pred[h] <- 0
              err[h] <- 0
            } else{
              sts_can <- sts[which(nams[h] == sts$ID), ]
              alts <- sts$ALT[match(dst, sts$ID)]
              lats <- sts$LAT[match(dst, sts$ID)]
              lons <- sts$LON[match(dst, sts$ID)]
              
              #probability of ocurrence prediction
              b <- as.numeric(ref)
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
              e <- sqrt(sum((yr - predict(fmt, type = 'response')) ^ 2)/(length(yr) - 3))
              e <- round((e * (maxi - mini)) + mini, 3)
              if(e < 0.001 & e > 0)
                e  <- 0.001
              
                #evaluating estimate
                wd_pred[h] <- pb
                raw_pred[h] <- p
                mod_pred[h] <- p
                err[h] <- e
                if(pb <= 0.5) mod_pred[h] <- 0
            }
          }
          if(intermediate){
          write.table(data.frame(ID = nams, 
                                 obs = as.numeric(x), 
                                 wd_pred = as.numeric(wd_pred), 
                                 raw_pred = as.numeric(raw_pred), 
                                 mod_pred = as.numeric(mod_pred),
                                 err = as.numeric(err)), 
                      output, quote =F, 
                      sep = '\t', row.names = F, na = '')}
          return(mod_pred)
        }
      }
  }
  
  
  # .fillData(x = prec[10,], datess, dat = 10, neibs = 10, distanc = distanc, sts = sts)
  