model.frame.censReg <- function( formula, ... ) {

   result <- model.frame.lm( formula )

   return( result )
}

