estfun.censReg <- function( x, includeSigma = TRUE, ... ) {

   result <- maxLik:::estfun.maxLik( x, ... )

   if( !includeSigma ) {
      result <- result[ , colnames( result ) != "logSigma" ]
   }

   return( result )
}

