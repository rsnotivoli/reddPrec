#Grid creation


gridPcp <- function (filled, points, sts, inidate, enddate, distanc){
  dir.create("./gridded/", showWarnings = F)
  
  datess = seq.Date(inidate, enddate, by = "day")
  predday = function(x, datess, filled, distanc = distanc, 
                     points = points, sts = sts) {
    d = filled[x, ]
    if (max(d, na.rm = T) == 0) {
      predpoint <- matrix(0, ncol = 2, nrow = nrow(points))
    } else {
      ff_pred <- function(h, distanc, d, sts, points, w){
        wsts <- match(names(distanc[[h]][[1]]), sts$ID)
        obs <- d[wsts]
        if (max(obs, na.rm = T) == 0) {
          return(c(0,0))
        } else {
          sts_can <- points[h, ]
          sts_nns <- sts[wsts, ]
          b <- obs
          b[b > 0] = 1
          y = b
          alt = sts_nns$ALT
          lat = sts_nns$Y
          lon = sts_nns$X
          fmtb <- suppressWarnings(glm(y ~ alt + lat + 
                                         lon, family = binomial()))
          ##
          newdata = data.frame(alt = sts_can$ALT, lat = sts_can$Y, 
                               lon = sts_can$X)
          pb <- predict(fmtb, newdata = newdata, type = "response")
          if (pb < 0.001 & pb > 0) pb <- 0.001
          pb <- round(pb, 3)
          mini = min(obs)/2
          maxi = max(obs) + (max(obs) - min(obs))
          yr = as.numeric((obs - mini)/(maxi - mini))
          fmt <- suppressWarnings(glm(yr ~ alt + lat + 
                                        lon, family = quasibinomial()))
          
          p <- predict(fmt, newdata = newdata, type = "response")
          if (p < 0.001 & p > 0) p <- 0.001
          p <- round((p * (maxi - mini)) + mini, 3)
          err <- sqrt(sum((yr - predict(fmt, type = "response"))^2)/(length(yr) - 
                                                                       3))
          if (err < 0.001 & err > 0) err <- 0.001
          err <- round((err * (maxi - mini)) + mini, 3)
          if (pb < 0.5) {
            ppred = 0
            eerr = 0
          } else if (pb >= 0.5 & p < 1) {
            ppred = 1
            eerr = err
          } else {
            ppred = p
            eerr = err
          }
        }
        return(c(ppred, eerr))
      }
      predpoint <- t(sfSapply(1:nrow(points), fun = ff_pred, 
                              distanc = distanc, d = d, sts = sts, 
                              points = points, w = w))
    }
    predpoint <- cbind(points$ID, predpoint)
    colnames(predpoint) = c("ID", "pred", "err")
    write.table(predpoint, paste("./gridded/", datess[x], ".txt", 
                                 sep = ""), quote = F, row.names = F, sep = "\t", 
                na = "")
  }
  for(dd in 1:length(datess)){
    print(paste('Computing date', datess[dd]))
    predday(x = dd, datess, filled, distanc = distanc, 
            points = points, sts = sts)
  }
}
