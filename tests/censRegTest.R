library( censReg )

data( "Affairs", package = "AER" )
affairsFormula <- affairs ~ age + yearsmarried + religiousness +
   occupation + rating

## usual tobit estimation
estResult <- censReg( affairsFormula, data = Affairs )
print.default( estResult )
print( estResult )
print( estResult, logSigma = FALSE )
maxLik:::summary.maxLik( estResult )
summary( estResult )
print( summary( estResult ), logSigma = FALSE )
coef( estResult )
coef( estResult, logSigma = FALSE )
vcov( estResult )
vcov( estResult, logSigma = FALSE )
coef( summary( estResult ) )
coef( summary( estResult ), logSigma = FALSE )
logLik( estResult )

## usual tobit estimation, BHHH method
estResultBhhh <- censReg( affairsFormula, data = Affairs, method = "BHHH" )
print.default( estResultBhhh )
print( estResultBhhh )
maxLik:::summary.maxLik( estResultBhhh )
summary( estResultBhhh )

## usual tobit estimation, BFGS method
estResultBfgs <- censReg( affairsFormula, data = Affairs, method = "BFGS" )
print.default( estResultBfgs )
print( estResultBfgs )
maxLik:::summary.maxLik( estResultBfgs )
summary( estResultBfgs )

## usual tobit estimation, NM method
estResultNm <- censReg( affairsFormula, data = Affairs, method = "NM" )
print.default( estResultNm )
print( estResultNm )
maxLik:::summary.maxLik( estResultNm )
summary( estResultNm )

## usual tobit estimation, SANN method
estResultSann <- censReg( affairsFormula, data = Affairs, method = "SANN" )
print.default( estResultSann )
print( estResultSann )
maxLik:::summary.maxLik( estResultSann )
summary( estResultSann )

## usual tobit estimation with user-defined starting values
estResultStart <- censReg( affairsFormula, data = Affairs,
   start = c( 8.17, -0.18, 0.55, -1.69, 0.33, -2.3, 2.13 ) )
print.default( estResultStart )
print( estResultStart )
maxLik:::summary.maxLik( estResultStart )
summary( estResultStart )
logLik( estResultStart )

## estimation with left-censoring at 5
Affairs$affairsAdd <- Affairs$affairs + 5
estResultAdd <- censReg( affairsAdd ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, left = 5 )
print.default( estResultAdd )
print( estResultAdd )
maxLik:::summary.maxLik( estResultAdd )
summary( estResultAdd )
coef( estResultAdd )
coef( estResultAdd, logSigma = FALSE )
vcov( estResultAdd )
vcov( estResultAdd, logSigma = FALSE )
logLik( estResultAdd )

## estimation with right-censoring
Affairs$affairsNeg <- - Affairs$affairs
estResultNeg <- censReg( affairsNeg ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, left = -Inf, right = 0 )
print.default( estResultNeg )
print( estResultNeg )
maxLik:::summary.maxLik( estResultNeg )
summary( estResultNeg )
coef( estResultNeg )
coef( estResultNeg, logSigma = FALSE )
vcov( estResultNeg )
vcov( estResultNeg, logSigma = FALSE )
logLik( estResultNeg )

## estimation with right-censoring at -5
Affairs$affairsAddNeg <- - Affairs$affairsAdd
estResultAddNeg <- censReg( affairsAddNeg ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, left = -Inf, right = -5 )
print.default( estResultAddNeg )
print( estResultAddNeg )
maxLik:::summary.maxLik( estResultAddNeg )
summary( estResultAddNeg )
coef( estResultAddNeg )
coef( estResultAddNeg, logSigma = FALSE )
vcov( estResultAddNeg )
vcov( estResultAddNeg, logSigma = FALSE )
logLik( estResultAddNeg )

## estimation with left and right censoring
estResultBoth <- censReg( affairsFormula, data = Affairs, right = 4 )
print.default( estResultBoth )
print( estResultBoth )
maxLik:::summary.maxLik( estResultBoth )
summary( estResultBoth )
print( summary( estResultBoth ), logSigma = FALSE )
coef( estResultBoth )
coef( estResultBoth, logSigma = FALSE )
vcov( estResultBoth )
vcov( estResultBoth, logSigma = FALSE )
coef( summary( estResultBoth ) )
coef( summary( estResultBoth ), logSigma = FALSE )
logLik( estResultBoth )
