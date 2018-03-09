
#' This function iterates the Quality Control by single time steps,
#' it can be parallelized. Apply after .qcFirst

#' @param x vector of precipitation (previously cleaned) of a single time step and an extra figure with the previous numebr of detected supects
#' @param y vector of precipitation (original) of a single time step and an extra figure with the previous numebr of detected supects
#' @param distanc matrix of distances with the names of the neighbouring stations
#' @param it integer number of iteration
#' @param sts data.frame with the stations info
#' @param datess dates vector in seq.Date format
#' @param printmeta Logical value wether metadata will be saved as a file

.qcLast <- function(x, y, dat, distanc, it, sts, datess, printmeta) {
  
  nams <- as.character(sts$ID) #names
  
  w_ori <- which(!is.na(y)) #available on original
  w_prec <- which(!is.na(x)) #available on cleaned
  if(length(w_ori) == 0 | length(w_prec) == 0) { #if no data, ends
    return(c(y, 0))
  } else{
    
    #subset dataset based on available data
    sub <- x[w_prec]
    or <- y
    dd <- distanc[match(nams, rownames(distanc)),]
    nam_d <- nams[w_prec] #stations with data in cleaned
    if (max(or, na.rm = T) == 0) {
      return(c(y, 0))
    }
      else {
        #start evaluating data
        code <- rep(NA, length(or))
        for (h in 1:length(or)) {
          can <- or[h]
          if(is.na(can))
              next else{
            #now using cleaned as neighbours
              #nearest stations with data
              wd <- which(nam_d == nams[h])
              if(length(wd) > 0){
                nr <- nam_d[-c(wd)]
                } else {nr <- nam_d}
              ww <- match(nr,dd[h,])
              ww <- sort(ww)[1:(neibs+1)]
              dst <- dd[h,][ww]
              ref <- sub[match(dst, nr)]
            
            #checks for codes
            if (max(ref, na.rm = T) == 0 & can > 0) {
              code[h] <- 1
            } else if (min(ref, na.rm = T) > 0 & can == 0) {
              code[h] <- 2
            } else{
              if(max(ref, na.rm = T) == 0){
                pb <- 0
                p <- 0
              }
              #compute glms
              else {
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
                mini <- min(ref, na.rm = T)/2
                maxi <- max(ref, na.rm = T) + (max(ref, na.rm = T) - min(ref, na.rm = T))
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
          }
        }
        
        #encoding detected
        ww = which(!is.na(code))
        if (length(ww) > 0) {
          k <- cbind(as.character(sts$ID)[ww], as.character(datess[dat]), 
                     code[ww], y[ww])
          colnames(k) = c("ID", "date", "code", "data")
          y[match(unique(nams[ww]), colnames(prec))] <- NA
          if (printmeta) {
            dir.create("./meta/", showWarnings = F)
            write.table(k, paste("./meta/meta_", as.character(datess[dat]), ".txt", 
                                 sep = ""), quote = F, row.names = F, na = "", 
                        sep = "\t")
          }
        }
        dim(y) <- c(output = length(y))
        return(y)
      }
  }
}
