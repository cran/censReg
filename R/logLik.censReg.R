logLik.censReg <- function( object, ... ) {

   result <- maxLik:::logLik.maxLik( object )
   attr( result, "df" ) <- sum( activePar( object ) )

   class( result ) <- "logLik"

   return( result )
}
