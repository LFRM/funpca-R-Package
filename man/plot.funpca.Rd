\name{plot.funpca}
\alias{plot.funpca}
\title{Plot fitted components
}
\description{
Plots fitted signals and shows acf/pacf for the each one. Additionally a plot for all curves is added at the beginning.
}
\usage{
\method{plot}{funpca}(x,...)
}
\arguments{
  \item{x}{\code{funpca} object.
}
  \item{\dots}{Other arguments to be called by plot(). 
}
}
\details{
Plot of the fitted results.
}
\value{
The function returns the selected plots. 
}
\references{
Rosales, F.\cr
work in progress
}
\author{
Francisco Rosales
}

\seealso{
\code{\link[funpca]{plot.funpca}} (package funpca)
}
\examples{
library(fda)
sdata <- NULL
data <- CanadianWeather$monthlyTemp
for(i in 1:ncol(data)) sdata <- cbind(sdata,spline(data[,i])$y)
x <- funpca(sdata, k=3)
plot(x)
}

\keyword{plot}

