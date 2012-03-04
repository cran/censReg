## log likelihood function for cross-sectional data
censRegLogLikCross <- function( beta, yVec, xMat, left, right, 
      obsBelow, obsBetween, obsAbove ) {
   yHat <- xMat %*% beta[ - length( beta ) ]
   sigma <- exp( beta[ length( beta ) ] )
   ll <- rep( NA, length( yVec ) )
   ll[ obsBelow ] <-
      pnorm( ( left - yHat[ obsBelow ] ) / sigma, log.p = TRUE )
   ll[ obsBetween ] <-
      dnorm( ( yVec - yHat )[ obsBetween ] / sigma, log = TRUE ) -
      log( sigma )
   ll[ obsAbove ] <-
      pnorm( ( yHat[ obsAbove ] - right ) / sigma, log.p = TRUE )

   ## gradients of log likelihood function for cross-sectional data
   grad <- matrix( NA, nrow = length( yVec ), ncol = length( beta ) )
   grad[ obsBelow, ] <-
      dnorm( ( left - yHat[ obsBelow ] ) / sigma ) /
      pnorm( ( left - yHat[ obsBelow ] ) / sigma ) *
      cbind( - xMat[ obsBelow, , drop = FALSE ] / sigma,
         - ( left - yHat[ obsBelow ] ) / sigma )
   grad[ obsBetween, ] <-
      cbind( ( ( yVec - yHat )[ obsBetween ] / sigma ) *
         xMat[ obsBetween, , drop = FALSE ] / sigma,
         ( ( yVec - yHat )[ obsBetween ] / sigma )^2 - 1 )
   grad[ obsAbove, ] <-
      dnorm( ( yHat[ obsAbove ] - right ) / sigma ) /
      pnorm( ( yHat[ obsAbove ] - right ) / sigma ) *
      cbind( xMat[ obsAbove, , drop = FALSE ] / sigma,
         - ( yHat[ obsAbove ] - right ) / sigma )
   attr( ll, "gradient" ) <- grad
   return( ll )
}

