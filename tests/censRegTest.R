library( censReg )
library( lmtest )
library( sandwich )

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
margEff( estResult )
summary( margEff( estResult ) )
logLik( estResult )
nobs( estResult )
extractAIC( estResult )
formula( estResult )
model.frame( estResult )
estfun( estResult )[ 20 * c(1:30), ]
meat( estResult )
bread( estResult )
sandwich( estResult )
all.equal( sandwich( estResult ), vcov( estResult ) )
waldtest( estResult, . ~ . - age )
waldtest( estResult, . ~ . - age, vcov = sandwich( estResult ) )

## usual tobit estimation, BHHH method
estResultBhhh <- censReg( affairsFormula, data = Affairs, method = "BHHH" )
print.default( estResultBhhh )
print( estResultBhhh )
margEff( estResultBhhh, returnJacobian = TRUE )
summary( margEff( estResultBhhh ) )
maxLik:::summary.maxLik( estResultBhhh )
summary( estResultBhhh )
all.equal( -crossprod( estfun( estResultBhhh ) ), 
   hessian( estResultBhhh ), check.attributes = FALSE )
all.equal( sandwich( estResultBhhh ), vcov( estResultBhhh ) )

## usual tobit estimation, BFGS method
estResultBfgs <- censReg( affairsFormula, data = Affairs, method = "BFGS" )
print.default( estResultBfgs )
print( estResultBfgs )
margEff( estResultBfgs, calcVCov = FALSE )
summary( margEff( estResultBfgs ) )
maxLik:::summary.maxLik( estResultBfgs )
summary( estResultBfgs )

## usual tobit estimation, NM method
estResultNm <- censReg( affairsFormula, data = Affairs, method = "NM" )
print.default( estResultNm )
print( estResultNm )
margEff( estResultNm )
summary( margEff( estResultNm ) )
maxLik:::summary.maxLik( estResultNm )
summary( estResultNm )

## usual tobit estimation, SANN method
estResultSann <- censReg( affairsFormula, data = Affairs, method = "SANN" )
print.default( estResultSann )
print( estResultSann )
margEff( estResultSann )
summary( margEff( estResultSann ) )
maxLik:::summary.maxLik( estResultSann )
summary( estResultSann )

## usual tobit estimation with user-defined starting values
estResultStart <- censReg( affairsFormula, data = Affairs,
   start = c( 8.17, -0.18, 0.55, -1.69, 0.33, -2.3, 2.13 ) )
print.default( estResultStart )
print( estResultStart )
margEff( estResultStart )
summary( margEff( estResultStart, calcVCov = FALSE, returnJacobian = TRUE ) )
maxLik:::summary.maxLik( estResultStart )
summary( estResultStart )
logLik( estResultStart )
nobs( estResultStart )
formula( estResultStart )

## estimation with left-censoring at 5
Affairs$affairsAdd <- Affairs$affairs + 5
estResultAdd <- censReg( affairsAdd ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, left = 5 )
print.default( estResultAdd )
print( estResultAdd )
margEff( estResultAdd )
summary( margEff( estResultAdd, returnJacobian = TRUE ) )
maxLik:::summary.maxLik( estResultAdd )
summary( estResultAdd )
coef( estResultAdd )
coef( estResultAdd, logSigma = FALSE )
vcov( estResultAdd )
vcov( estResultAdd, logSigma = FALSE )
logLik( estResultAdd )
nobs( estResultAdd )
extractAIC( estResultAdd )

## estimation with right-censoring
Affairs$affairsNeg <- - Affairs$affairs
estResultNeg <- censReg( affairsNeg ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, left = -Inf, right = 0 )
print.default( estResultNeg )
print( estResultNeg )
margEff( estResultNeg, calcVCov = FALSE, returnJacobian = TRUE )
summary( margEff( estResultNeg ) )
maxLik:::summary.maxLik( estResultNeg )
summary( estResultNeg )
coef( estResultNeg )
coef( estResultNeg, logSigma = FALSE )
vcov( estResultNeg )
vcov( estResultNeg, logSigma = FALSE )
logLik( estResultNeg )
nobs( estResultNeg )
extractAIC( estResultNeg )
model.frame( estResultNeg )

## estimation with right-censoring at -5
Affairs$affairsAddNeg <- - Affairs$affairsAdd
estResultAddNeg <- censReg( affairsAddNeg ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, left = -Inf, right = -5 )
print.default( estResultAddNeg )
print( estResultAddNeg )
margEff( estResultAddNeg )
summary( margEff( estResultAddNeg, calcVCov = FALSE ) )
maxLik:::summary.maxLik( estResultAddNeg )
summary( estResultAddNeg )
coef( estResultAddNeg )
coef( estResultAddNeg, logSigma = FALSE )
vcov( estResultAddNeg )
vcov( estResultAddNeg, logSigma = FALSE )
logLik( estResultAddNeg )
nobs( estResultAddNeg )
extractAIC( estResultAddNeg )

## estimation with left and right censoring
estResultBoth <- censReg( affairsFormula, data = Affairs, right = 4 )
print.default( estResultBoth )
print( estResultBoth )
margEff( estResultBoth )
summary( margEff( estResultBoth ) )
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
nobs( estResultBoth )
extractAIC( estResultBoth )
estfun( estResultBoth )[ 20 * c(1:30), ]
meat( estResultBoth )
bread( estResultBoth )
sandwich( estResultBoth )
all.equal( sandwich( estResultBoth ), vcov( estResultBoth ) )
waldtest( estResultBoth, . ~ . - age )
waldtest( estResultBoth, . ~ . - age, vcov = sandwich( estResultBoth ) )

## with empty levels
Affairs2 <- Affairs
Affairs2$religiousness <- as.factor( Affairs2$religiousness )
Affairs2 <- Affairs2[ Affairs2$religiousness != "5", ]
estResultEmpty <- censReg( affairsFormula, data = Affairs2 )
print.default( estResultEmpty )
print( estResultEmpty )
summary( estResultEmpty )
coef( estResultEmpty )
vcov( estResultEmpty )
margEff( estResultEmpty )
summary( margEff( estResultEmpty ) )
formula( estResultEmpty )
model.frame( estResultEmpty )
estfun( estResultEmpty )[ 20 * c(1:26), ]
meat( estResultEmpty )
bread( estResultEmpty )
sandwich( estResultEmpty )
all.equal( sandwich( estResultEmpty ), vcov( estResultEmpty ) )
waldtest( estResultEmpty, . ~ . - age )
waldtest( estResultEmpty, . ~ . - age, vcov = sandwich( estResultEmpty ) )


# returning log-likelihood contributions only (no estimations)
logLikBhhh <- censReg( affairsFormula, data = Affairs, method = "BHHH",
   start = coef( estResultBhhh ), logLikOnly = TRUE )
print( logLikBhhh )
all.equal( sum( logLikBhhh ), c( logLik( estResultBhhh ) ) )
logLikStart <- censReg( affairsFormula, data = Affairs,
   start = c( 8.17, -0.18, 0.55, -1.69, 0.33, -2.3, 2.13 ),
   logLikOnly = TRUE )
print( logLikStart )

