\name{coef.summary.censReg}
\alias{coef.summary.censReg}

\title{Coefficients of Censored Regression Models and their Statistical Properties}

\description{
   This function returns the estimated coefficients 
   of censored regression models
   as well as their standard errors, z-values, and P-values.   
}

\usage{
\method{coef}{summary.censReg}( object, logSigma = TRUE, ... )
}

\arguments{
   \item{object}{object of class \code{"summary.censReg"}
      (returned by \code{\link{summary.censReg}}).}
   \item{logSigma}{logical value indicating whether the variance(s)
      of the model should be returned logarithmized.}
   \item{\dots}{currently not used.}
}
\value{
   \code{coef.summary.censReg} returns an matrix,
   where each row corresponds to one coefficient
   and the 4 rows contain the estimated coefficients,
   their standard errors, z-values, and P-values.
}

\author{Arne Henningsen}

\seealso{\code{\link{censReg}}, \code{\link{summary.censReg}} and 
   \code{\link{coef.censReg}}}

\examples{
## Kleiber & Zeileis ( 2008 ), page 142
data( "Affairs", package = "AER" )
estResult <- censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs )
coef( summary( estResult ) )
}

\keyword{methods}
