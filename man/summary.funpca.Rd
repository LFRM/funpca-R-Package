\name{summary.funpca}
\alias{summary.funpca}
\title{funpca Summary}
\description{Takes an \code{funpca} object produced by \code{funpca} and summarizes the information of the components fi (individual fits). 
}
\usage{
\method{summary}{funpca}(x,...)
}
\arguments{
  \item{x}{\code{funpca} object. 
  }
  \item{\dots}{further arguments to be passed to summary().
} 
}
\value{
The function gives basic statistics of the components resulting from applying \code{funpca}. 
}
\references{
Rosales, F. and Krivobokova, T.\cr
work in progress. 
}
\author{
Francisco Rosales
}

\seealso{
\code{\link[funpca]{plot.funpca}} (package funpca),
}
\examples{
library(fda)
sdata <- NULL
data <- CanadianWeather$monthlyTemp
for(i in 1:ncol(data)) sdata <- cbind(sdata,spline(data[,i])$y)
x <- funpca(sdata, k=3)
}

\keyword{summary}
\keyword{fda}
