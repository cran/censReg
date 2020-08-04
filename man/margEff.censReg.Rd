\name{margEff.censReg}
\alias{margEff.censReg}
\alias{summary.margEff.censReg}

\title{Marginal Effects in Censored Regression Models}

\description{
   The \code{margEff} method computes the marginal effects of the explanatory variables
   on the expected value of the dependent variable evaluated.
   Please note that this functionality is currently not available
   for panel data models.
}

\usage{
\method{margEff}{censReg}( object, xValues = NULL, vcov = NULL,
   calcVCov = TRUE, returnJacobian = FALSE, ... )

\method{summary}{margEff.censReg}( object, ... )
}

\arguments{
   \item{object}{argument \code{object} of the \code{margEff} method
      must be an object of class \code{"censReg"}
      (returned by \code{\link{censReg}});
      argument \code{object} of the \code{summary} method
      must be an object of class \code{"margEff.censReg"}
      (returned by \code{margEff.censReg}).}
   \item{xValues}{vector that specifies the values of the explanatory variables
      (including the intercept if it is included in the model),
      at which the marginal effects should be calculated.
      The number and order of the elements of this vector
      must correspond to the number and order of the estimated coefficients
      (without sigma).
      If this argument is \code{NULL} (or not specified),
      argument \code{xValues} is set to the mean values
      of the explanatory variables.}
   \item{vcov}{a symmetric matrix that specifies the variance covariance
      matrix of the estimated coefficients (without \code{logSigma})
      that should be used to calculate the variance covariance matrix
      and the standard errors of the marginal effects.
      If this argument is \code{NULL} (the default),
      the variance covariance matix is obtained by
      \code{vcov( object, logSigma = FALSE )}.} 
   \item{calcVCov}{logical. If \code{TRUE}, 
     the approximate variance covariance matrices of the marginal effects 
     is calculated and returned as an attribute (see below).}
   \item{returnJacobian}{logical. If \code{TRUE}, 
      the Jacobian of the marginal effects with respect to the coefficients 
      is returned as an attribute (see below).}
   \item{\dots}{currently not used.}
}

\value{
   \code{margEff.censReg} returns an object of class \code{"margEff.censReg"},
   which is a vector of the marginal effects of the explanatory variables
   on the expected value of the dependent variable evaluated
   at the mean values of the explanatory variables.
   The returned object has an attribute \code{df.residual}, 
   which is equal to the degrees of freedom of the residuals. 
   
   If argument \code{calcVCov} is \code{TRUE}, 
   the object returned by \code{margEff.censReg} has an attribute \code{vcov}, 
   which is a the approximate variance covariance matrices 
   of the marginal effects calculated
   with the Delta method.
   
   If argument \code{returnJacobian} is \code{TRUE}, 
   the object returned by \code{margEff.censReg}has an attribute \code{jacobian}, 
   which is the Jacobian of the marginal effects 
   with respect to the coefficients.

   \code{summary.margEff.censReg} returns 
   an object of class \code{"summary.margEff.censReg"},
   which is a matrix with four columns,
   where the first column contains the marginal effects,
   the second column contains the standard errors of the marginal effects,
   the third column contains the corresponding t-values,
   and the fourth columns contains the corresponding P values.
}

\author{Arne Henningsen}

\seealso{\code{\link{censReg}}, \code{\link{coef.censReg}}, 
   and \code{\link{summary.censReg}}}

\examples{
## Kleiber & Zeileis ( 2008 ), page 142
data( "Affairs", package = "AER" )
estResult <- censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs )
margEff( estResult )
summary( margEff( estResult ) )
margEff( estResult, xValues = c( 1, 40, 4, 2, 4, 4 ) )
}

\keyword{methods}
