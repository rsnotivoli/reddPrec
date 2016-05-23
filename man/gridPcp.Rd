\name{gridPcp}
\alias{gridPcp}

\title{
Gridding daily precipitation.
}
\description{
This function uses serially-completed data to estimate new predicted values in each point and day of a given point dataset.
}
\usage{
gridPcp(filled, points, sts, inidate, enddate, parallel = TRUE, ncpu = 2)
}

\arguments{
  \item{filled}{
Object of class \code{matrix} containing the daily precipitation data. This dataset can be the \code{Filled.RData} resulting from \code{gapFilling} function. Each column represents one station. The names of columns have to be names of the stations.
}
  \item{points}{
Object of class \code{matrix} containing the points of the new grid info. Must have at least four fields: \code{ID}: station identifier; \code{ALT}: altitude; \code{X}: Longitude in UTM projection (meters); and \code{Y}: Latitude in UTM projection (meters). Tabulation separated.
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
Logical. If \code{TRUE}, parallel computing is activated and the processes will be distributed among the \code{ncpu} number of processor cores.
}
  \item{ncpu}{
Only if \code{parallel = TRUE}. Sets the number of processor cores used to parallel computing.
}
}
\details{
The precipitation data used to make the grid should be a serially-completed dataset without any missing value. If we use a dataset with gaps, the number of near stations to each point will change along the period and we will introduce inhomogeneities to the final series. We recommend to use the output of \code{gapFilling} function.
}
\value{
A new \code{gridded} directory will be created in working path with one file per day. Each file contains a \code{data.frame} with many rows as points. The columns show the identifier (\code{ID})of each station, the magnitude prediction \code{pred} and the standard error \code{err} of the model. All of these values are expresed in the same units as original values (mm*10).
}
\author{
Roberto Serrano-Notivoli
}
\note{
The time computing depends directly on the number of points of the \code{points} dataset. If you need to fill a large grid we recommend to use a powerful computer with many number of proccesor cores (to distribute properly the work) and RAM (to process the whole amount of data).
}

\examples{

  #loads example data
  data(precipDataset)

  #runs function (only for two days)
  gridPcp(filled=filled,points=points,sts=sts,inidate=as.Date('2001-01-01'),
    enddate=as.Date('2001-01-01'),parallel=TRUE,ncpu=2)

}
