\name{summary.funpca}
\alias{summary.funpca}
\title{funpca Summary}
\description{Takes an \code{funpca} object produced by \code{funpca} and summarizes the information of the components fi (individual fits). 
}
\usage{
\method{summary}{funpca}(object,...)
}
\arguments{
  \item{object}{\code{funpca} object. 
  }
  \item{\dots}{further arguments to be passed to summary().
} 
}
\value{
The function gives basic statistics of the components resulting from applying \code{funpca}. 
}
\references{
Rosales, F. and Krivobokova, T.\cr
For more details see <https://ediss.uni-goettingen.de/handle/11858/00-1735-0000-0028-87F9-6>  
}
\author{
Francisco Rosales \email{francisco.rosales-marticorena@protonmail.com}
}

\seealso{
\code{\link[funpca]{plot.funpca}} (package funpca),
}
\examples{
\donttest{
library(fda)
sdata <- NULL
data <- CanadianWeather$monthlyTemp
for(i in 1:ncol(data)) sdata <- cbind(sdata,spline(data[,i])$y)
x <- funpca(sdata, k=3)
summary(x)
}
}

\keyword{summary}
\keyword{fda}
