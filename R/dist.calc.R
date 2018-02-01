#Distance calculation
  #returns a list of n elements (same length as the 
  #number of points grid) with the 10 nearest station
  #names in each one.


dist.calc <- function(filled, sts, points, ncpus){
  x1 <- cbind(sts$X, sts$Y)
  x2 <- cbind(points$X, points$Y)
  distanc <- rdist(x1, x2)/1000
  colnames(distanc) <- points$ID
  rownames(distanc) <- sts$ID
  #check complete filled dataset
  wna <- length(which(colSums(is.na(filled)) > 0))
  
  if(wna > 0) warning('The input dataset contains NAs, selected 
                      stations in some days could contain NAs')
    ff_10nearest <- function(x, nams){
      kk <- data.frame(N = nams, D = x)
      kk <- kk[order(kk$D),]
      kk <- kk[1:10,]
      dists <- kk[,2]
      names(dists) <- kk[,1]
      list(dists)
    }
    print('Computing distances')
    sfInit(parallel = TRUE, cpus = ncpus)
    distanc <- sfApply(distanc, 2, ff_10nearest, nams = rownames(distanc))
    
  distanc
  save(distanc, file = 'distances.RData')
}