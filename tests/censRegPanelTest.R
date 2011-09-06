library( censReg )
library( plm )

nId <- 15
nTime <- 4

set.seed( 123 )
pData <- data.frame(
   id = rep( paste( "F", 1:nId, sep = "_" ), each = nTime ),
   time = rep( 1980 + 1:nTime, nId ) )
pData$ui <- rep( rnorm( nId ), each = nTime )
pData$x1 <- rnorm( nId * nTime )
pData$x2 <- runif( nId * nTime )
pData$ys <- -1 + pData$ui + 2 * pData$x1 + 3 * pData$x2 + rnorm( nId * nTime )
pData$y <- ifelse( pData$ys > 0, pData$ys, 0 )
nData <- pData # save data set without information on panel structure
pData <- pdata.frame( pData, c( "id", "time" ) )


## Newton-Raphson method
randEff <- censReg( y ~ x1 + x2, data = pData )
print( randEff )
print( randEff, logSigma = FALSE )
maxLik:::summary.maxLik( randEff )
summary( randEff )
print( summary( randEff ), logSigma = FALSE )
coef( randEff )
coef( randEff, logSigma = FALSE )
vcov( randEff )
vcov( randEff, logSigma = FALSE )
coef( summary( randEff ) )
coef( summary( randEff ), logSigma = FALSE )
logLik( randEff )
extractAIC( randEff )
print.default( randEff )


## BHHH method
randEffBhhh <- censReg( y ~ x1 + x2, data = pData, method = "BHHH" )
print( randEffBhhh )
maxLik:::summary.maxLik( randEffBhhh )
summary( randEffBhhh )
print.default( randEffBhhh )


## BFGS method (optim)
randEffBfgs <- censReg( y ~ x1 + x2, data = pData, method = "BFGS" )
print( randEffBfgs )
maxLik:::summary.maxLik( randEffBfgs )
summary( randEffBfgs )
print.default( randEffBfgs )


## BFGS method (R)
randEffBfgsr <- censReg( y ~ x1 + x2, data = pData, method = "BFGSR" )
print( randEffBfgsr )
maxLik:::summary.maxLik( randEffBfgsr )
summary( randEffBfgsr )
print.default( randEffBfgsr )


## BHHH with starting values
randEffBhhhStart <- censReg( y ~ x1 + x2, data = pData, method = "BHHH",
   start = c( -0.4, 1.7, 2.2, -0.1, -0.01 ) )
print( randEffBhhhStart )
summary( randEffBhhhStart )


## left-censoring at 5
pData$yAdd <- pData$y + 5
randEffAdd <- censReg( yAdd ~ x1 + x2, data = pData, method = "BFGSR", left = 5 )
print( randEffAdd )
maxLik:::summary.maxLik( randEffAdd )
summary( randEffAdd )
coef( randEffAdd )
coef( randEffAdd, logSigma = FALSE )
vcov( randEffAdd )
vcov( randEffAdd, logSigma = FALSE )
logLik( randEffAdd )
extractAIC( randEffAdd )
print.default( randEffAdd )


## right-censoring
pData$yNeg <- - pData$y
randEffNeg <- censReg( yNeg ~ x1 + x2, data = pData, method = "BFGSR",
   left = -Inf, right = 0 )
print( randEffNeg )
maxLik:::summary.maxLik( randEffNeg )
summary( randEffNeg )
coef( randEffNeg )
coef( randEffNeg, logSigma = FALSE )
vcov( randEffNeg )
vcov( randEffNeg, logSigma = FALSE )
logLik( randEffNeg )
extractAIC( randEffNeg )
print.default( randEffNeg )


## right-censoring at -5
pData$yAddNeg <- - pData$yAdd
randEffAddNeg <- censReg( yAddNeg ~ x1 + x2, data = pData, method = "BFGSR",
   left = -Inf, right = -5 )
print( randEffAddNeg )
maxLik:::summary.maxLik( randEffAddNeg )
summary( randEffAddNeg )
coef( randEffAddNeg )
coef( randEffAddNeg, logSigma = FALSE )
vcov( randEffAddNeg )
vcov( randEffAddNeg, logSigma = FALSE )
logLik( randEffAddNeg )
extractAIC( randEffAddNeg )
print.default( randEffAddNeg )


## both right and left censoring
pData$yBoth <- ifelse( pData$y < 3, pData$y, 3 )
randEffBoth <- censReg( yBoth ~ x1 + x2, data = pData, method = "BFGSR",
   left = 0, right = 3 )
print( randEffBoth )
maxLik:::summary.maxLik( randEffBoth )
summary( randEffBoth )
print( summary( randEffBoth ), logSigma = FALSE )
coef( randEffBoth )
coef( randEffBoth, logSigma = FALSE )
vcov( randEffBoth )
vcov( randEffBoth, logSigma = FALSE )
coef( summary( randEffBoth ) )
coef( summary( randEffBoth ), logSigma = FALSE )
logLik( randEffBoth )
extractAIC( randEffBoth )
print.default( randEffBoth )


## re-order observations/individuals
set.seed( 234 )
perm <- sample( nId )
nData2 <- nData
nData2$id <- NA
for( i in 1:nId ) {
   nData2$id[ nData$id == paste( "F", i, sep = "_" ) ] <-
      paste( "G", perm[ i ], sep = "_" )
}
pData2 <- pdata.frame( nData2, c( "id", "time" ) )
randEffBfgsr2 <- censReg( y ~ x1 + x2, data = pData2, method = "BFGSR" )
all.equal( randEffBfgsr2[ -c(11,12) ], randEffBfgsr[ -c(11,12) ] )
all.equal( sort( randEffBfgsr2[[ 11 ]] ), sort( randEffBfgsr[[ 11 ]] ) )


## unbalanced panel data
nDataUnb <- nData[ -c( 2, 5, 6, 8 ), ]
pDataUnb <- pdata.frame( nDataUnb, c( "id", "time" ) )
randEffBfgsrUnb <- censReg( y ~ x1 + x2, data = pDataUnb, method = "BFGSR" )
print( randEffBfgsrUnb )
maxLik:::summary.maxLik( randEffBfgsrUnb )
summary( randEffBfgsrUnb )
logLik( randEffBfgsrUnb )
extractAIC( randEffBfgsrUnb )
print.default( randEffBfgsrUnb )


## NAs in data
pDataNa <- pData
obsNa <- which( ! rownames( pData ) %in% rownames( pDataUnb ) )
pDataNa$y[ obsNa[ 1:2 ] ] <- NA
pDataNa$x1[ obsNa[ 3 ] ] <- NA
pDataNa$x2[ obsNa[ c( 1, 2, 4 ) ] ] <- NA
randEffBfgsrNa <- censReg( y ~ x1 + x2, data = pDataNa, method = "BFGSR" )
all.equal( randEffBfgsrNa[ -12 ], randEffBfgsrUnb[ -12 ] )


