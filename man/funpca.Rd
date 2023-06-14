\name{funpca}
\alias{funpca}
\title{
Functional Principal Component Analysis
}
\description{
Performs functional principal component analysis using the mixed models representation of smoothing splines. 
}
\usage{
funpca(mat,k)
}
\arguments{
  \item{mat}{Is a rectangular matrix with no missing values.Each colum represents a sample.  
}
  \item{k}{Desired number of eigen functions to construct subj spec deviations. Should be between 1 and the sample size.  
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
\item{residuals}{Remainder component for each subject.}
\item{y}{Data used for all the computations.}
\item{call}{Call of funpca.}
}
\references{
Rosales, F. \cr
For more details see <https://ediss.uni-goettingen.de/handle/11858/00-1735-0000-0028-87F9-6> 
}
\author{
Francisco Rosales \email{francisco.rosales-marticorena@protonmail.com},
}
\seealso{
\code{\link[fda]{fda}} (package fda)
}
\examples{
\donttest{
library(fda)
sdata <- NULL
data <- CanadianWeather$monthlyTemp
for(i in 1:ncol(data)) sdata <- cbind(sdata,spline(data[,i])$y)
x <- funpca(sdata, k=3)
}
}

\keyword{nonlinear}
\keyword{funpca}
\keyword{splines}