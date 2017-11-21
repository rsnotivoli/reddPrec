\name{points}
\alias{points}
\docType{data}
\title{
Information about points of the grid.
}
\description{
Information of the all points related to \code{precip} dataset with an identifier, longitude, latitude and altitude.
}
\usage{data("precipDataset")}
\format{
  A \code{data.frame} with 676 points on the following 4 variables.
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
