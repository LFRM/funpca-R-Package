\name{funpca-package}
\alias{funpca-package}
\docType{package}
\title{
Functional Principal Component Analysis
}
\description{
Performs functional principal component analysis using the mixed models representation of smoothing splines. 
}
\details{
\tabular{ll}{
Package: \tab fpcamm\cr
Version: \tab 1.0\cr
Date: \tab 2023-06-08\cr
Depends: \tab Brobdingnag, MASS, nlme\cr
}

Index:
\preformatted{
funpca                    Performs FPCA using the MM representation of
                          penalized splines.
plot.funpca   	          Plots fitted curves: overall trend, subj spec
                          deviations and derivative of the mean curve.
summary.funpca            Summary of individual fits.
}
The function \code{funpca()} is used to fit the model. Using the resulting \code{funpca} object, fitted curves or their derivatives can be plotted with \code{\link{plot}} and summary information on the fit can be printed using \code{\link{summary}}.
}

\author{
Francisco Rosales
Maintainer: Francisco Rosales \email{francisco.rosales-marticorena@protonmail.com}
}

\references{
Rosales, F.\cr
For more details see <https://ediss.uni-goettingen.de/handle/11858/00-1735-0000-0028-87F9-6> 
}

\keyword{package}
\seealso{
\code{\link[fda]{fda}} (package fda)
}
