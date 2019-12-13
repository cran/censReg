library( "censReg" )
library( "lmtest" )
library( "sandwich" )

options( digits = 5 )

printAll <- function( x, logSigmaFalse = FALSE, sDigits = 2, 
      meCalcVCov = TRUE, meReturnJacobian = FALSE,
      sumMeCalcVCov = TRUE, sumMeReturnJacobian = FALSE ) {
   for( n in names( x ) ) {
      if( ! n %in% c( "code", "message", "iterations" ) ) {
         cat( "$", n, "\n", sep = "" )
         if( n %in% c( "estimate", "hessian", "gradientObs" ) ) {
            print( round( x[[ n ]], 2 ) )
         } else if( n %in% c( "gradient" ) ) {
            print( round( x[[ n ]], 3 ) )
         } else {
            print( x[[ n ]] )
         }
         cat( "\n" )
      }
   }
   cat( "class\n" )
   print( class( x ) )
   
   cat( "print( x, digits = 2 )\n" )
   print( x, digits = 2 )
   
   if( logSigmaFalse ) {
      cat( "print( x, logSigma = FALSE, digits = 2 )\n" )
      print( x, logSigma = FALSE, digits = 2 )
   }
   
   cat( "print( round( margEff( x ), digits = 2 ) )\n" )
   print( round( margEff( x, calcVCov = meCalcVCov, 
      returnJacobian = meReturnJacobian ), 2 ) )
   
   cat( "printME( margEff( x ) )\n" )
   printME( margEff( x, calcVCov = meCalcVCov,
      returnJacobian = meReturnJacobian ) )
   
   cat( "print( summary( margEff( x ) ), digits = sDigits )\n" )
   print( summary( margEff( x, calcVCov = sumMeCalcVCov, 
      returnJacobian = sumMeReturnJacobian ) ), digits = sDigits )
   
   x$code <- 0
   x$message <- "removed message"
   x$iterations <- 0
   
   cat( "print( maxLik:::summary.maxLik( x ), sDigits )\n" )
   print( maxLik:::summary.maxLik( x ), digits = sDigits )
   
   cat( "print( summary( x ), digits = sDigits )\n" )
   print( summary( x ), digits = sDigits )
   
   if( logSigmaFalse ) {
      cat( "print( summary( x ), logSigma = FALSE, digits = sDigits )\n" )
      print( summary( x ), logSigma = FALSE, digits = sDigits )
   }
}

printME <- function( x ) {
   print( round( x, 3 ) )
   y <- attributes( x )
   for( n in names( y ) ) {
      if( ! n %in% c( "names" ) ) {
         cat( "attr(,\"", n, "\")\n", sep = "" )
         if( n %in% c( "vcov" ) ) {
            print( round( y[[ n ]], 3 ) )
         } else if( n %in% c( "jacobian" ) ) {
            print( round( y[[ n ]], 3 ) )
         } else {
            print( y[[ n ]] )
         }
      }
   }
}

data( "Affairs", package = "AER" )
affairsFormula <- affairs ~ age + yearsmarried + religiousness +
   occupation + rating

## usual tobit estimation
estResult <- censReg( affairsFormula, data = Affairs )
printAll( estResult, logSigmaFalse = TRUE, sDigits = 1 )
round( coef( estResult ), 2 )
round( coef( estResult, logSigma = FALSE ), 2 )
round( vcov( estResult ), 2 )
round( vcov( estResult, logSigma = FALSE ), 2 )
round( coef( summary( estResult ) ), 2 )
round( coef( summary( estResult ), logSigma = FALSE ), 2 )
all.equal( margEff( estResult ),
   margEff( estResult, xValues = estResult$xMean ) )
round( margEff( estResult, xValues = c( 1, 40, 4, 2, 4, 4 ) ), 2 )
printME( margEff( estResult, xValues = c( 1, 40, 4, 2, 4, 4 ) ) )
print( summary( margEff( estResult, xValues = c( 1, 40, 4, 2, 4, 4 ) ) ),
   digits = 2 )
logLik( estResult )
nobs( estResult )
extractAIC( estResult )
formula( estResult )
model.frame( estResult )
round( estfun( estResult )[ 20 * c(1:30), ], 2 )
round( meat( estResult ), 2 )
round( bread( estResult ), 2 )
round( sandwich( estResult ), 2 )
# all.equal( sandwich( estResult ), vcov( estResult ) )
waldtest( estResult, . ~ . - age )
waldtest( estResult, . ~ . - age, vcov = sandwich( estResult ) )

## usual tobit estimation, BHHH method
estResultBhhh <- censReg( affairsFormula, data = Affairs, method = "BHHH" )
printAll( estResultBhhh, meReturnJacobian = TRUE )
all.equal( -crossprod( estfun( estResultBhhh ) ), 
   hessian( estResultBhhh ), check.attributes = FALSE )
all.equal( sandwich( estResultBhhh ), vcov( estResultBhhh ) )

## usual tobit estimation, BFGS method
estResultBfgs <- censReg( affairsFormula, data = Affairs, method = "BFGS" )
printAll( estResultBfgs, meCalcVCov = FALSE )

## usual tobit estimation, NM method
estResultNm <- censReg( affairsFormula, data = Affairs, method = "NM" )
printAll( estResultNm )

## usual tobit estimation, SANN method
estResultSann <- censReg( affairsFormula, data = Affairs, method = "SANN" )
printAll( estResultSann )

## usual tobit estimation with user-defined starting values
estResultStart <- censReg( affairsFormula, data = Affairs,
   start = c( 8.17, -0.18, 0.55, -1.69, 0.33, -2.3, 2.13 ) )
printAll( estResultStart, sumMeCalcVCov = FALSE, sumMeReturnJacobian = TRUE )
logLik( estResultStart )
nobs( estResultStart )
formula( estResultStart )

## estimation with left-censoring at 5
Affairs$affairsAdd <- Affairs$affairs + 5
estResultAdd <- censReg( affairsAdd ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, left = 5 )
printAll( estResultAdd, sumMeReturnJacobian = TRUE )
round( coef( estResultAdd ), 2 )
round( coef( estResultAdd, logSigma = FALSE ), 2 )
round( vcov( estResultAdd ), 2 )
round( vcov( estResultAdd, logSigma = FALSE ), 2 )
logLik( estResultAdd )
nobs( estResultAdd )
extractAIC( estResultAdd )

## estimation with right-censoring
Affairs$affairsNeg <- - Affairs$affairs
estResultNeg <- censReg( affairsNeg ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, left = -Inf, right = 0 )
printAll( estResultNeg, meCalcVCov = FALSE, meReturnJacobian = TRUE )
round( coef( estResultNeg ), 2 )
round( coef( estResultNeg, logSigma = FALSE ), 2 )
round( vcov( estResultNeg ), 2 )
round( vcov( estResultNeg, logSigma = FALSE ), 2 )
logLik( estResultNeg )
nobs( estResultNeg )
extractAIC( estResultNeg )
model.frame( estResultNeg )

## estimation with right-censoring at -5
Affairs$affairsAddNeg <- - Affairs$affairsAdd
estResultAddNeg <- censReg( affairsAddNeg ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, left = -Inf, right = -5 )
printAll( estResultAddNeg, sumMeCalcVCov = FALSE )
round( coef( estResultAddNeg ), 2 )
round( coef( estResultAddNeg, logSigma = FALSE ), 2 )
round( vcov( estResultAddNeg ), 2 )
round( vcov( estResultAddNeg, logSigma = FALSE ), 2 )
logLik( estResultAddNeg )
nobs( estResultAddNeg )
extractAIC( estResultAddNeg )

## estimation with left and right censoring
estResultBoth <- censReg( affairsFormula, data = Affairs, right = 4 )
printAll( estResultBoth, logSigmaFalse = TRUE )
round( coef( estResultBoth ), 2 )
round( coef( estResultBoth, logSigma = FALSE ), 2 )
round( vcov( estResultBoth ), 2 )
round( vcov( estResultBoth, logSigma = FALSE ), 2 )
round( coef( summary( estResultBoth ) ), 2 )
round( coef( summary( estResultBoth ), logSigma = FALSE ), 2 )
logLik( estResultBoth )
nobs( estResultBoth )
extractAIC( estResultBoth )
round( estfun( estResultBoth )[ 20 * c(1:30), ], 2 )
round( meat( estResultBoth ), 2 )
round( bread( estResultBoth ), 2 )
round( sandwich( estResultBoth ), 2 )
# all.equal( sandwich( estResultBoth ), vcov( estResultBoth ) )
waldtest( estResultBoth, . ~ . - age )
waldtest( estResultBoth, . ~ . - age, vcov = sandwich( estResultBoth ) )

## with empty levels
Affairs2 <- Affairs
Affairs2$religiousness <- as.factor( Affairs2$religiousness )
Affairs2 <- Affairs2[ Affairs2$religiousness != "5", ]
estResultEmpty <- censReg( affairsFormula, data = Affairs2 )
printAll( estResultEmpty )
round( coef( estResultEmpty ), 2 )
round( vcov( estResultEmpty ), 2 )
formula( estResultEmpty )
model.frame( estResultEmpty )
round( estfun( estResultEmpty )[ 20 * c(1:26), ], 2 )
round( meat( estResultEmpty ), 2 )
round( bread( estResultEmpty ), 1 )
round( sandwich( estResultEmpty ), 2 )
# all.equal( sandwich( estResultEmpty ), vcov( estResultEmpty ) )
waldtest( estResultEmpty, . ~ . - age )
waldtest( estResultEmpty, . ~ . - age, vcov = sandwich( estResultEmpty ) )


# returning log-likelihood contributions only (no estimations)
logLikBhhh <- censReg( affairsFormula, data = Affairs, method = "BHHH",
   start = coef( estResultBhhh ), logLikOnly = TRUE )
round( c( logLikBhhh ), 2 )
round( attr( logLikBhhh, "gradient" ), 2 )
all.equal( sum( logLikBhhh ), c( logLik( estResultBhhh ) ) )
logLikStart <- censReg( affairsFormula, data = Affairs,
   start = c( 8.17, -0.18, 0.55, -1.69, 0.33, -2.3, 2.13 ),
   logLikOnly = TRUE )
round( c( logLikStart ), 2 )
round( attr( logLikStart, "gradient" ), 2 )

