\name{qcPrec}
\alias{qcPrec}
\title{
Quality control: Identifies and removes suspect data from the original dataset
}
\description{
This function uses original data to estimate new predicted values and compare them with observations. If exist large differences, removes the original values.
}
\usage{
qcPrec(prec, sts, inidate, enddate, parallel = TRUE,
      ncpu = 2, printmeta = TRUE, thres = NA)
}
\arguments{
  \item{prec}{
Object of class \code{matrix} containing the original precipitation data. Each column represents one station. The names of columns have to be names of the stations.
}
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
Logical. When \code{TRUE}, parallel computing is activated and the processes will be distributed among the \code{ncpu} number of processor cores.
}
  \item{ncpu}{
Only when \code{parallel = TRUE}. Sets the number of processor cores used to parallel computing.
}
  \item{printmeta}{
When \code{TRUE}, one file per day will be written in subdirectory \code{./meta}.
}
\item{thres}{
Threshold applied to search nearest stations. If \code{thres=NA} the function will search 10 nearest observations without a distance limit. A positive number indicates the threshold in kilometers.
}
}
\details{
The process of quality control uses five criteria to flag suspect data. All of them are based on the calculation of reference values (RV) made with the 10 nearest observations (NNS) that day. For this reason, a minimum of 11 available data by day is mandatory. The five criteria are : 1) Suspect data: Observed > 0 & all their 10 NNS == 0; 2) Suspect zero: Observed == 0 & all their 10 NNS > 0; 3) Suspect outlier: Observed is 10 times higher or lower than RV; 4) Suspect wet: Observed == 0, wet probability is over 99\%, and predicted magnitude is over 5 litres and 5) Suspect dry: Observed > 5 litres, dry probability is over 99\%, and predicted magnitude is under 0.1 litres.

All of these criteria are prepared to work with precipitation in tenths (milimetres*10).
}
\value{
A new file called \code{cleaned.RData} will be created in working directory. The load of this file (\code{load('cleaned.RData')}) will add a matrix with the original data filtered by quality control. If \code{printmeta = TRUE}, a new \code{meta} directory will be created in working path with one file per day. Each file contains a \code{data.frame} with many rows as flagged data in that day. The columns show the identifier (\code{ID})of each station; the \code{date}; the criteria \code{code} through the data was flagged and the removed \code{data}. There are five different codes referred to the five criteria: 1 = Suspect data; 2 = Suspect zero; 3 = Suspect outlier; 4 = Suspect wet and 5 = Suspect dry.
}
\author{
Roberto Serrano-Notivoli
}
\note{
It is possible that many \code{warnings} appear at the end of the execution of the function. Probably most of them will be due to a bad fitting in algorithm used. This is because some days in some locations, the nearest stations are so different that they are not enough to explain strongly the precipitation in one location.
}
\examples{
  #loads example data
  data(precipDataset)

  #runs function
  qcPrec(prec=preci,sts=sts,inidate=as.Date('2001-01-01'),
    enddate=as.Date('2001-01-02'),parallel=TRUE,ncpu=2,printmeta=TRUE,thres=NA)
}
