summary.censReg <- function( object, ... ) {

   result <- maxLik:::summary.maxLik( object )

   result$call <- object$call
   result$nObs <- object$nObs

   class( result ) <- c( "summary.censReg", class( result ) )

   return( result )
}
