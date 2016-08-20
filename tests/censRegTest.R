library( "censReg" )
library( "lmtest" )
library( "sandwich" )

options( digits = 5 )

printAll <- function( x ) {
   for( n in names( x ) ) {
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
   cat( "class\n" )
   print( class( x ) )
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
printAll( estResult )
print( estResult, digits = 2 )
print( estResult, logSigma = FALSE, digits = 2 )
print( maxLik:::summary.maxLik( estResult ), digits = 1 )
print( summary( estResult ), digits = 1 )
print( summary( estResult ), logSigma = FALSE, digits = 1 )
round( coef( estResult ), 2 )
round( coef( estResult, logSigma = FALSE ), 2 )
round( vcov( estResult ), 2 )
round( vcov( estResult, logSigma = FALSE ), 2 )
round( coef( summary( estResult ) ), 2 )
round( coef( summary( estResult ), logSigma = FALSE ), 2 )
round( margEff( estResult ), 2 )
printME( margEff( estResult ) )
summary( margEff( estResult ) )
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
meat( estResult )
round( bread( estResult ), 2 )
round( sandwich( estResult ), 2 )
# all.equal( sandwich( estResult ), vcov( estResult ) )
waldtest( estResult, . ~ . - age )
waldtest( estResult, . ~ . - age, vcov = sandwich( estResult ) )

## usual tobit estimation, BHHH method
estResultBhhh <- censReg( affairsFormula, data = Affairs, method = "BHHH" )
printAll( estResultBhhh )
print( estResultBhhh, digits = 2 )
round( margEff( estResultBhhh, returnJacobian = TRUE ), 2 )
printME( margEff( estResultBhhh, returnJacobian = TRUE ) )
print( summary( margEff( estResultBhhh ) ), digits = 2 )
print( maxLik:::summary.maxLik( estResultBhhh ), digits = 2 )
print( summary( estResultBhhh ), digits = 2 )
all.equal( -crossprod( estfun( estResultBhhh ) ), 
   hessian( estResultBhhh ), check.attributes = FALSE )
all.equal( sandwich( estResultBhhh ), vcov( estResultBhhh ) )

## usual tobit estimation, BFGS method
estResultBfgs <- censReg( affairsFormula, data = Affairs, method = "BFGS" )
printAll( estResultBfgs )
print( estResultBfgs, digits = 2 )
round( margEff( estResultBfgs, calcVCov = FALSE ), 2 )
printME( margEff( estResultBfgs, calcVCov = FALSE ) )
print( summary( margEff( estResultBfgs ) ), digits = 2 )
print( maxLik:::summary.maxLik( estResultBfgs ), 2 )
print( summary( estResultBfgs ), digits = 2 )

## usual tobit estimation, NM method
estResultNm <- censReg( affairsFormula, data = Affairs, method = "NM" )
printAll( estResultNm )
print( estResultNm, digits = 2 )
round( margEff( estResultNm ), 2 )
printME( margEff( estResultNm ) )
print( summary( margEff( estResultNm ) ), digits = 2 )
print( maxLik:::summary.maxLik( estResultNm ), digits = 2 )
print( summary( estResultNm ), digits = 2 )

## usual tobit estimation, SANN method
estResultSann <- censReg( affairsFormula, data = Affairs, method = "SANN" )
printAll( estResultSann )
print( estResultSann, digits = 2 )
round( margEff( estResultSann ), 2 )
printME( margEff( estResultSann ) )
print( summary( margEff( estResultSann ) ), digits = 2 )
print( maxLik:::summary.maxLik( estResultSann ), digits = 2 )
print( summary( estResultSann ), digits = 2 )

## usual tobit estimation with user-defined starting values
estResultStart <- censReg( affairsFormula, data = Affairs,
   start = c( 8.17, -0.18, 0.55, -1.69, 0.33, -2.3, 2.13 ) )
printAll( estResultStart )
print( estResultStart, digits = 2 )
round( margEff( estResultStart ), 2 )
printME( margEff( estResultStart ) )
print( summary(
   margEff( estResultStart, calcVCov = FALSE, returnJacobian = TRUE ) ),
   digits = 2 )
print( maxLik:::summary.maxLik( estResultStart ), digits = 2 )
print( summary( estResultStart ), digits = 2 )
logLik( estResultStart )
nobs( estResultStart )
formula( estResultStart )

## estimation with left-censoring at 5
Affairs$affairsAdd <- Affairs$affairs + 5
estResultAdd <- censReg( affairsAdd ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, left = 5 )
printAll( estResultAdd )
print( estResultAdd, digits = 2 )
round( margEff( estResultAdd ), 2 )
printME( margEff( estResultAdd ) )
print( summary( margEff( estResultAdd, returnJacobian = TRUE ) ), digits = 2 )
print( maxLik:::summary.maxLik( estResultAdd ), digits = 2 )
print( summary( estResultAdd ), digits = 2 )
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
printAll( estResultNeg )
print( estResultNeg, digits = 2 )
round( margEff( estResultNeg, calcVCov = FALSE, returnJacobian = TRUE ), 2 )
printME( margEff( estResultNeg, calcVCov = FALSE, returnJacobian = TRUE ) )
print( summary( margEff( estResultNeg ) ), digits = 2 )
print( maxLik:::summary.maxLik( estResultNeg ), digits = 2 )
print( summary( estResultNeg ), digits = 2 )
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
printAll( estResultAddNeg )
print( estResultAddNeg, digits = 2 )
round( margEff( estResultAddNeg ), 2 )
printME( margEff( estResultAddNeg ) )
print( summary( margEff( estResultAddNeg, calcVCov = FALSE ) ), digits = 2 )
print( maxLik:::summary.maxLik( estResultAddNeg ), digits = 2 )
print( summary( estResultAddNeg ), digits = 2 )
round( coef( estResultAddNeg ), 2 )
round( coef( estResultAddNeg, logSigma = FALSE ), 2 )
round( vcov( estResultAddNeg ), 2 )
round( vcov( estResultAddNeg, logSigma = FALSE ), 2 )
logLik( estResultAddNeg )
nobs( estResultAddNeg )
extractAIC( estResultAddNeg )

## estimation with left and right censoring
estResultBoth <- censReg( affairsFormula, data = Affairs, right = 4 )
printAll( estResultBoth )
print( estResultBoth, digits = 2 )
round( margEff( estResultBoth ), 2 )
printME( margEff( estResultBoth ) )
print( summary( margEff( estResultBoth ) ), digits = 2 )
print( maxLik:::summary.maxLik( estResultBoth ), digits = 2 )
print( summary( estResultBoth ), digits = 2 )
print( summary( estResultBoth ), logSigma = FALSE, digits = 2 )
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
print( estResultEmpty, digits = 2 )
print( summary( estResultEmpty ), digits = 2 )
round( coef( estResultEmpty ), 2 )
round( vcov( estResultEmpty ), 2 )
round( margEff( estResultEmpty ), 2 )
printME( margEff( estResultEmpty ) )
print( summary( margEff( estResultEmpty ) ), 2 )
formula( estResultEmpty )
model.frame( estResultEmpty )
round( estfun( estResultEmpty )[ 20 * c(1:26), ], 2 )
round( meat( estResultEmpty ), 2 )
round( bread( estResultEmpty ), 2 )
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

