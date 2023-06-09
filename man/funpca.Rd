\name{funpca}
\alias{funpca}
\title{
Functional Principal Component Analysis
}
\description{
Performs functional principal component analysis using the mixed models representation of smoothing splines. 
}
\usage{
funpca(DATA,k,correlation=NULL)
}
\arguments{
  \item{DATA}{Is a rectangular matrix with no missing values.Each colum represents a sample.  
}
  \item{k}{Desired number of eigen functions to construct subj spec deviations. Should be between 1 and the sample size.  
}
  \item{correlation}{Is an entry of the form: corARMA(p=P,q=Q). It represents the correlation to be used for the remainder component. If left empty a white noise remainder is assumed. 
}
}
\details{
The method assumes DATA is a complete rectangular matrix and hence does not support missing values. 
}
\value{
A list object of class \code{funpca} containing the following information.
\item{est}{Mixed model estimation}
\item{f}{A matrix with the fitted overall trend. All columns contain the same information}
\item{di}{A matrix with the fitted subj spec deviations}
\item{fi}{Fitted values for each subject, i.e. fitted overall trend + fitted subj spec deviations + subj spec seasonality.}
\item{error}{Remainder component for each subject.}
\item{fd1}{x first derivative of the overall trend.}
\item{fd2}{x second derivative of the overall trend.}
\item{residuals}{Remainder component for each subject.}
\item{DATA}{Data used for all the computations.}
\item{call}{Call of funpca.}
}
\references{
Rosales, F. \cr
work in progress
}
\author{
Francisco Rosales \email{l.rosalesmarticore at stud.uni-goettingen.de},
}
\seealso{
\code{\link[fda]{fda}} (package fda)
}
\examples{
library(fda)
sdata <- NULL
data <- CanadianWeather$monthlyTemp
for(i in 1:ncol(data)) sdata <- cbind(sdata,spline(data[,i])$y)
x <- funpca(sdata, k=3)
}

\keyword{nonlinear}
\keyword{funpca}
\keyword{splines}