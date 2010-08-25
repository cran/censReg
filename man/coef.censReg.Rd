\name{coef.censReg}
\alias{coef.censReg}
\alias{vcov.censReg}
\alias{logLik.censReg}

\title{Coefficients, their Covariances, and Log-Likelihood Values of Censored Regression Models}

\description{
   These functions extract the coefficient vectors,
   the corresponding covariance matrices,
   and log-likelihood values
   from censored regression models.
}

\usage{
\method{coef}{censReg}( object, logSigma = TRUE, ... )

\method{vcov}{censReg}( object, logSigma = TRUE, ... )

\method{logLik}{censReg}( object, ... )
}

\arguments{
   \item{object}{object of class \code{"censReg"}
      (returned by \code{\link{censReg}}).}
   \item{logSigma}{logical value indicating whether the variance(s)
      of the model should be returned logarithmized.}
   \item{\dots}{currently not used.}
}
\value{
   \code{coef.censReg} returns a vector of the estimated coefficients.

   \code{vcov.censReg} returns the covariance matrix
   of the estimated coefficients.

   \code{logLik.censReg} returns an object of class \code{"logLik"}.
   This object is the log-likelihood value of the estimated model
   and has an attribute \code{"df"}
   that gives the degrees of freedom,
   i.e.\ the number of estimated parameters.
}

\note{
   When the censored regression model is estimated,
   the log-likelihood function is maximized with respect
   to the coefficients and the \emph{logarithm(s)}
   of the variance(s).
   Hence, if argument \code{logSigma} is \code{FALSE},
   the variance(s) of the model is/are calculated by applying
   the exponential function to the estimated logarithmized variance(s)
   and the covariance matrix of all parameters is calculated
   by the Delta method.
}

\author{Arne Henningsen}

\seealso{\code{\link{censReg}}, \code{\link{summary.censReg}}, 
   and \code{\link{coef.summary.censReg}}}

\examples{
## Kleiber & Zeileis ( 2008 ), page 142
data( "Affairs", package = "AER" )
estResult <- censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs )
coef( estResult )
coef( estResult, logSigma = FALSE )
vcov( estResult )
vcov( estResult, logSigma = FALSE )
logLik( estResult )
}

\keyword{methods}
