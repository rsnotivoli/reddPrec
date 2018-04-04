
#' This function Uses the great-circle distances (Haversine formula),
#' returns a vector with the ID of the sorted neighbours

#' @param x vector with the coordinates of the candidate locations c(lat,lon) 
#' @param y matrix or data.frame with the coordinates of the reference locations (including candidate or not) and their IDs two columns [lat,lon,ID]
#' @param thres maximum radius where neighbouring stations will be searched

.dist_near <- function(x, y, thres){
  
  #convert degrees in radians
  lat <- as.numeric(x[1])*pi/180
  lon <- as.numeric(x[2])*pi/180
  latN <- as.numeric(y[,1])*pi/180
  lonN <- as.numeric(y[,2])*pi/180
  nams <- as.character(y[,3])
  #compute distances
  distN <- sqrt((sin((lat-latN)/2))^2 +
                  cos(lat)*cos(latN)*(sin((lon-lonN)/2))^2)
  Rx2 <- 2*6371 # earth mean radius (km)
  distN <- (Rx2) * asin(pmin(distN,1))
  w <- which(distN == 0)
  
  if(length(w) > 0){
    #select only the 10 nearest points
    #and return their indices
    names(distN) <- 1:length(distN)
    if (!is.na(thres)){
      distN <- sort(distN)
      distN <- distN[distN<thres]
      distN <- as.numeric(names(distN)[2:length(distN)])
      distN <- nams[distN]
    }
    else{
      distN <- sort(distN)
      distN <- as.numeric(names(distN)[2:length(distN)])
      distN <- nams[distN]
    }
  } else{
    names(distN) <- 1:length(distN)
    if (!is.na(thres)){
      distN <- sort(distN)
      distN <- distN[distN<thres]
      distN <- as.numeric(names(distN))
      distN <- nams[distN]
    }
    else{
      distN <- sort(distN)
      distN <- as.numeric(names(distN))
      distN <- nams[distN]
    }
    }
  
  return(distN)
}
