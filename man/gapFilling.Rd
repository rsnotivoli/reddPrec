\name{gapFilling}
\alias{gapFilling}
\title{
Daily precipitation reconstruction. Makes predictions to all days and locations.
}
\description{
This function uses original data to estimate new predicted values in each location and day.
}
\usage{
gapFilling(prec, sts, inidate, enddate, parallel = TRUE,
          ncpu = 2, thres = NA)
}
\arguments{
  \item{prec}{
Object of class \code{matrix} containing the daily precipitation data. This dataset can be the \code{cleaned.RData} resulting from \code{qcPrec} function. Each column represents one station. The names of columns have to be names of the stations.}
  \item{sts}{
Object of class \code{matrix} containing the stations info. Must have at least four fields: \code{ID}: station identifier; \code{ALT}: altitude; \code{X}: Longitude in UTM projection (meters); and \code{Y}: Latitude in UTM projection (meters). Tabulation separated.
}
  \item{inidate}{
Object of class \code{Date} in format \code{'YYYY-mm-dd'} defining the first day of quality control process
}
  \item{enddate}{
Object of class \code{Date} in format \code{'YYYY-mm-dd'} defining the last day of quality control process
}
  \item{parallel}{
Logical. If \code{TRUE}, parallel computing is activated and the processes will be distributed among the \code{ncpu} number of processor cores.
}
  \item{ncpu}{
Only if \code{parallel = TRUE}. Sets the number of processor cores used to parallel computing.
}
\item{thres}{
Threshold applied to search nearest stations. If \code{thres=NA} the function will search 10 nearest observations without a distance limit. A positive number indicates the threshold in kilometers.
}
}
\value{
A new file called \code{Filled.RData} will be created in working directory. The load of this file (\code{load('Filled.RData')}) will add a matrix containing the original data with missing values filled with this function. A new \code{days} directory will be created in working path with one file per day. Each file contains a \code{data.frame} with many rows as stations. The columns show the identifier (\code{ID})of each station; the observed value \code{obs} (if exists); the binomial prediction \code{predb} (dry<0.5; wet>=0.5); the raw prediction of precipitation magnitude \code{pred1} resulting from models; the corrected magnitude prediction \code{pred2} by binomial prediction (if dry, pred2==0); the final standardized magnitude prediction \code{pred3} and the standard error \code{err} of the model. All variables except \code{predb} are in the same units as original values (mm*10).
}
\author{
Roberto Serrano-Notivoli
}
\examples{

  #loads example data
  data(precipDataset)

  #runs function
  gapFilling(prec=cleaned,sts=sts,inidate=as.Date('2001-01-01'),
    enddate=as.Date('2001-01-02'),parallel=TRUE,ncpu=2,thres=NA)
}
