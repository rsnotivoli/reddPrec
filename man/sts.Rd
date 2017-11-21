\name{sts}
\alias{sts}
\docType{data}
\title{
Information about stations
}
\description{
Information of the 48 stations related to \code{precip} dataset with an identifier, longitude, latitude and altitude.
}
\usage{data("precipDataset")}
\format{
  A \code{data.frame} with 48 observations on the following 4 variables.
  \describe{
    \item{\code{ID}}{identifier of the station.}
    \item{\code{X}}{longitude of the station.}
    \item{\code{Y}}{latitude of the station.}
    \item{\code{ALT}}{altitude of the station.}
  }
}
\examples{
data(precipDataset)
}
\keyword{datasets}
