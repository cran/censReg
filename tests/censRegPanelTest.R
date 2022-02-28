library( "censReg" )
library( "plm" )

# load outputs that were previously produced by this script 
saved <- new.env()
load( "censRegPanelTest.RData.save", envir = saved )

options( digits = 5 )

printAll <- function( objName, what = "print" ) {
   cat( "Comparing new object '", objName, "' to previously saved object...",
      sep = "" )
   x <- get( objName )
   if( !exists( objName, envir = saved, inherits = FALSE ) ) {
      cat( " previously saved object not found\n" )
   } else {
      xSaved <- get( objName, envir = saved, inherits = FALSE )
      if( !isTRUE( all.equal( class( x ), class( xSaved ) ) ) ) {
         cat( " different classes:\n" )
         cat( "new:\n" )
         print( class( x ) )
         cat( "saved:\n" )
         print( class( xSaved ) )
      } else if( !isTRUE( all.equal( names( x ), names( xSaved ) ) ) ) {
         cat( " different names:\n" )
         cat( "new:\n" )
         print( names( x ) )
         cat( "saved:\n" )
         print( names( xSaved ) )
      } else {
         cat( "\n" )
      }
      for( n in names( x ) ) {
         if( ! n %in% c( "code", "gradient", "iterations", "last.step",
               "message" ) ) {
            cat( "   comparing component '", n, "' ...", sep = "" )
            if( n == "vcov" ) {
               tol <- 5e-1
            } else if( n == "estimate" ) {
               tol <- 5e-2
            } else {
               tol <- 5e-3
            }
            testRes <-  all.equal( x[[ n ]], xSaved[[ n ]], tol = tol )
            if( isTRUE( testRes ) ) {
               cat( " OK\n" )
            } else {
               cat( " different\n" )
               print( testRes )
               cat( "new:\n" )
               print( x[[ n ]] )
               cat( "saved:\n" )
               try( print( xSaved[[ n ]] ) )
            }
         }
      }
   }
   
   for( mName in c( "Coef", "CoefNoLs", "Vcov", "VcovNoLs",
         "CoefSum", "CoefSumNoLs", "LogLik", "Nobs", "ExtractAIC" ) ) {
      cat( "   comparing method '", mName, "' ...", sep = "" )
      tol <- 5e-3
      if( mName == "Coef" ) {
         xm <- coef( x )
         tol <- 5e-2
      } else if( mName == "CoefNoLs" ) {
         xm <- coef( x, logSigma = FALSE )
         tol <- 5e-2
      } else if( mName == "Vcov" ) {
         xm <- vcov( x )
         tol <- 5e-1
      } else if( mName == "VcovNoLs" ) {
         xm <- vcov( x, logSigma = FALSE )
         tol <- 5e-1
      } else if( mName == "CoefSum" ) {
         xm <- coef( summary( x ) )
         tol <- 5e-2
      } else if( mName == "CoefSumNoLs" ) {
         xm <- coef( summary( x ), logSigma = FALSE )
         tol <- 5e-2
      } else if( mName == "LogLik" ) {
         xm <- logLik( x )
      } else if( mName == "Nobs" ) {
         xm <- nobs( x )
      } else if( mName == "ExtractAIC" ) {
         xm <- extractAIC( x )
      } else {
         stop( "unknown value of 'mName': ", mName )
      }
      methodObjName <- paste0( objName, mName )
      if( !exists( methodObjName, envir = saved, inherits = FALSE ) ) {
         cat( " previously saved object not found\n" )
      } else {
         xmSaved <- get( methodObjName, envir = saved, inherits = FALSE )
         testRes <- all.equal( xm, xmSaved, tol = tol )
         if( isTRUE( testRes ) ) {
            cat( " OK\n" )
         } else {
            cat( " different\n" )
            print( testRes )
            cat( "new:\n" )
            print( xm )
            cat( "saved:\n" )
            print( xmSaved )
         }
      }
      # assign to parent frame so that it will be included in the saved workspace
      assign( methodObjName, xm, envir = parent.frame() )
   }
      
   if( what %in% c( "print", "methods", "all" ) ) {
      print( x, digits = 1 )
      print( x, logSigma = FALSE , digits = 1 )
      print( maxLik:::summary.maxLik( x ), digits = 1 )
      print( summary( x ), digits = 1 )
      print( summary( x ), logSigma = FALSE , digits = 1 )
   }
   if( what %in% c( "methods", "all" ) ) {
      print( round( coef( x ), 2 ) )
      print( round( coef( x, logSigma = FALSE ), 2 ) )
      print( round( vcov( x ), 2 ) )
      print( round( vcov( x, logSigma = FALSE ), 2 ) )
      print( round( coef( summary( x ) ), 2 ) )
      print( round( coef( summary( x ), logSigma = FALSE ), 2 ) )
      try( margEff( x ) )
      print( logLik( x ) )
      print( nobs( x ) )
      print( extractAIC( x ) )
   }
   
   if( what == "all" ) {
      for( n in names( x ) ) {
         cat( "$", n, "\n", sep = "" )
         if( n %in% c( "estimate", "gradientObs" ) ) {
            print( round( x[[ n ]], 2 ) )
         } else if( n %in% c( "hessian" ) ) {
            print( round( x[[ n ]], 1 ) )
         } else if( n %in% c( "gradient" ) ) {
         } else if( ! n %in% c( "last.step" ) ) {
            print( x[[ n ]] )
         }
         cat( "\n" )
      }
      cat( "class\n" )
      print( class( x ) )
   }
}

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
pData <- pdata.frame( pData, c( "id", "time" ), stringsAsFactors = FALSE )


## Newton-Raphson method
randEff <- censReg( y ~ x1 + x2, data = pData )
printAll( "randEff" )
try( margEff( randEff ) )
# only intercept
randEffOnlyInt <- censReg( y ~ 1, data = pData )
printAll( "randEffOnlyInt", what = "diff" )
# no intercept
randEffNoInt <- censReg( y ~ x1 -1, data = pData )
printAll( "randEffNoInt" )
# neither intercept nor explanatory variables
try( censReg( y ~ -1, data = pData ) )

## BHHH method
randEffBhhh <- censReg( y ~ x1 + x2, data = pData, method = "BHHH" )
printAll( "randEffBhhh" )


## BFGS method (optim)
randEffBfgs <- censReg( y ~ x1 + x2, data = pData, method = "BFGS" )
printAll( "randEffBfgs" )


## BFGS method (R)
randEffBfgsr <- censReg( y ~ x1 + x2, data = pData, method = "BFGSR" )
printAll( "randEffBfgsr", what = "none" )


## BHHH with starting values
randEffBhhhStart <- censReg( y ~ x1 + x2, data = pData, method = "BHHH",
   start = c( -0.4, 1.7, 2.2, -0.1, -0.01 ) )
printAll( "randEffBhhhStart" )


## left-censoring at 5
pData$yAdd <- pData$y + 5
randEffAdd <- censReg( yAdd ~ x1 + x2, data = pData, left = 5 )
printAll( "randEffAdd" )


## right-censoring
pData$yNeg <- - pData$y
randEffNeg <- censReg( yNeg ~ x1 + x2, data = pData, 
   left = -Inf, right = 0 )
printAll( "randEffNeg" )


## right-censoring at -5
pData$yAddNeg <- - pData$yAdd
randEffAddNeg <- censReg( yAddNeg ~ x1 + x2, data = pData, 
   left = -Inf, right = -5 )
printAll( "randEffAddNeg" )


## both right and left censoring
pData$yBoth <- ifelse( pData$y < 3, pData$y, 3 )
randEffBoth <- censReg( yBoth ~ x1 + x2, data = pData, 
   left = 0, right = 3 )
printAll( "randEffBoth" )


## re-order observations/individuals
set.seed( 234 )
perm <- sample( nId )
nData2 <- nData
nData2$id <- NA
for( i in 1:nId ) {
   nData2$id[ nData$id == paste( "F", i, sep = "_" ) ] <-
      paste( "G", perm[ i ], sep = "_" )
}
pData2 <- pdata.frame( nData2, c( "id", "time" ), stringsAsFactors = FALSE )
randEff2 <- censReg( y ~ x1 + x2, data = pData2 )
all.equal( randEff2[ -c(3,5,6,7,9,11,15) ],
   randEff[ -c(3,5,6,7,9,11,15) ], tolerance = 1e-2 )

# check if the order of observations/individuals influences the likelihood values
d1c1 <- censReg( y ~ x1 + x2, data = pData, start = coef(randEff),
   iterlim = 0 )
all.equal( d1c1[-c(5,6,7,9,12,15,19)], randEff[-c(5,6,7,9,12,15,19)] )
round( d1c1$maximum -  randEff$maximum, 12 )

d2c2 <- censReg( y ~ x1 + x2, data = pData2, start = coef(randEff2),
   iterlim = 0 )
all.equal( d2c2[-c(5,6,7,9,12,15,19)], randEff2[-c(5,6,7,9,12,15,19)] )
round( d2c2$maximum -  randEff2$maximum, 12 )

d1c2 <- censReg( y ~ x1 + x2, data = pData,  
   start = coef(randEff2), iterlim = 0 )
round( d2c2$maximum - d1c2$maximum, 12 )
round( d2c2$gradient - d1c2$gradient, 12 )

d2c1 <- censReg( y ~ x1 + x2, data = pData2, 
   start = coef(randEff), iterlim = 0 )
round( d1c1$maximum - d2c1$maximum, 12 )
round( d1c1$gradient - d2c1$gradient, 12 )

round( d2c2$maximum - d2c1$maximum, 3 )
round( d1c1$maximum - d1c2$maximum, 3 )

d1cS <- censReg( y ~ x1 + x2, data = pData, 
   start = randEff$start, iterlim = 0 )
d2cS <- censReg( y ~ x1 + x2, data = pData2, 
   start = randEff$start, iterlim = 0 )
round( d1cS$maximum - d2cS$maximum, 12 )
round( d1cS$gradient - d2cS$gradient, 12 )


## unbalanced panel data
nDataUnb <- nData[ -c( 2, 5, 6, 8 ), ]
pDataUnb <- pdata.frame( nDataUnb, c( "id", "time" ), stringsAsFactors = FALSE )
randEffUnb <- censReg( y ~ x1 + x2, data = pDataUnb )
printAll( "randEffUnb" )


## NAs in data
pDataNa <- pData
obsNa <- which( ! rownames( pData ) %in% rownames( pDataUnb ) )
pDataNa$y[ obsNa[ 1:2 ] ] <- NA
pDataNa$x1[ obsNa[ 3 ] ] <- NA
pDataNa$x2[ obsNa[ c( 1, 2, 4 ) ] ] <- NA
randEffNa <- censReg( y ~ x1 + x2, data = pDataNa )
all.equal( randEffNa[ -15 ], randEffUnb[ -15 ] )


# returning log-likelihood contributions only (no estimations)
logLikRandEff <- censReg( y ~ x1 + x2, data = pData, start = coef( randEff ),
   logLikOnly = TRUE )
print( logLikRandEff, digits = 1 )
all.equal( sum( logLikRandEff ), c( logLik( randEff ) ) )
logLikStart <- censReg( y ~ x1 + x2, data = pData, 
   start = c( -0.4, 1.7, 2.2, -0.1, -0.01 ), logLikOnly = TRUE )
print( round( c( logLikStart ), 3 ) )
print( round( attr( logLikStart, "gradient" ), 2 ) )


# save all objectives that were produced in this script
# (in order to compare them with objects created by this script in the future)
rm( saved )
save.image( "censRegPanelTest.RData" )

