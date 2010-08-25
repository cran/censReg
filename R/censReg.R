censReg <- function( formula, left = 0, right = Inf,
      data = sys.frame( sys.parent() ), start = NULL,
      nGHQ = 8, ... ) {

   ## checking formula
   if( class( formula ) != "formula" ) {
      stop( "argument 'formula' must be a formula" )
   } else if( length( formula ) != 3 ) {
      stop( "argument 'formula' must be a 2-sided formula" )
   }

   ## checking limits
   # left
   if( !is.numeric( left ) ) {
      stop( "argument 'left' must be a number" )
   } else if( length( left ) != 1 ) {
      stop( "argument 'left' must be a scalar (single number)" )
   }
   # right
   if( !is.numeric( right ) ) {
      stop( "argument 'right' must be a number" )
   } else if( length( right ) != 1 ) {
      stop( "argument 'right' must be a scalar (single number)" )
   }
   # both
   if( left >= right ) {
      stop( "argument 'right' must be a larger number than argument 'left'" )
   }

   ## preparing model matrix and model response
   mc <- match.call( expand.dots = FALSE )
   m <- match( "data", names( mc ), 0 )
   mf <- mc[ c( 1, m ) ]
   mf$formula <- formula
   attributes( mf$formula ) <- NULL
   mf$na.action <- na.pass
   mf[[ 1 ]] <- as.name( "model.frame" )
   mf <- eval( mf, parent.frame() )
   mt <- attr( mf, "terms" )
   xMat <- model.matrix( mt, mf )
   xNames <- colnames( xMat )
   yVec <- model.response( mf )
   yName <- as.character( formula )[ 2 ]
   if( length( yVec ) != nrow( xMat ) ) {
      stop( "the number of observations of the endogenous variable (",
         length( yVec ), ") is not equal to the number of observations",
         " of the exogenous variables (", nrow( xMat ), ")" )
   }

   ## extract information on panel structure of data set
   isPanel <- "pdata.frame" %in% class( data )
   if( isPanel ) {
      pIndex <- attributes( data )$index
      ## check if observations are ordered with respect to names of individuals
      # (theoretically, it is not required that the observations are ordered
      # alphabetically with respect to individuals' names but
      # if the individuals are not in the same order for each time period,
      # the observations are allocated to a wrong individual)
      if( !identical( order( pIndex[[ 1 ]] ), 1:length( pIndex[[ 1 ]] ) ) ) {
         stop( "names of individuals in attributes(data)$index[[1]]",
            " must be in alphabetical order but they are not;",
            " please fix this and re-run censReg()." )
      }

   }

   ## check if endogenous variable is within limits
   if( any( yVec[ !is.na( yVec ) ] < left ) ) {
      warning( "at least one value of the endogenous variable is smaller than",
         " the left limit" )
   } else if( any( yVec[ !is.na( yVec ) ] > right ) ) {
      warning( "at least one value of the endogenous variable is larger than",
         " the right limit" )
   }

   ## detect and remove observations with NAs, NaNs, and INFs
   validObs <- rowSums( is.na( cbind( yVec, xMat ) ) |
      is.infinite( cbind( yVec, xMat ) ) ) == 0
   yVec <- yVec[ validObs ]
   xMat <- xMat[ validObs, , drop = FALSE ]
   if( isPanel ) {
      pIndex <- pIndex[ validObs, , drop = FALSE ]
      indNames <- unique( pIndex[[ 1 ]] )  # 'names' of individuals
      nInd <- length( indNames )           # number of individuals
      timeNames <- unique( pIndex[[ 2 ]] ) # 'names' of time periods
      nTime <- length( timeNames )         # number of time periods
   }

   ## starting values
   nParam <- ncol( xMat ) + 1
   if( is.null( start ) ) {
      if( isPanel ) {
         assign( "validObs2", validObs, inherits = TRUE )
         # Random effects panel model estimation for starting values
         rEff <- plm( formula, data = data, subset = validObs2,
            effect = "individual", model = "random" )
         start <- c( coef( rEff ),
            0.5 * log( rEff$ercomp$sigma$id ),
            0.5 * log( rEff$ercomp$sigma$idios ) )
      } else {
         # OLS estimation for starting values
         ols <- lm.fit( xMat, yVec )
         start <- c( ols$coefficients,
            log( sum( ols$residuals^2 ) / length( ols$residuals ) ) )
      }
   } else {
      if( !is.numeric( start ) ) {
         stop( "argument 'start' must be numeric" )
      } else if( length( start ) != nParam ) {
         stop( "argument 'start' must have length ", nParam )
      }
   }

   if( isPanel ) {
      ## naming coefficients
      names( start ) <- c( colnames( xMat ), "logSigmaMu", "logSigmaNu" )

      ## Abscissae and weights for the Gauss-Hermite-Quadrature
      ghqPoints <- ghq( nGHQ, modified = FALSE )

      ## re-organize data 
      xArr <- array( NA, dim = c( nInd, nTime, ncol( xMat ) ) )
      yMat <- matrix( NA, nrow = nInd, ncol = nTime )
      for( i in 1:nTime ) {
         obsTime <- pIndex[[ 2 ]] == timeNames[ i ]
         xArr[ indNames %in% pIndex[[ 1 ]][ obsTime ], i, ] <- xMat[ obsTime, ]
         yMat[ indNames %in% pIndex[[ 1 ]][ obsTime ], i ] <- yVec[ obsTime ]
      }

      ## classify observations
      obsBelow <- yMat <= left & !is.na( yMat )
      obsAbove <- yMat >= right & !is.na( yMat )
      obsBetween <- !obsBelow & !obsAbove & !is.na( yMat )

      ## log likelihood function for panel data (incl. gradients)
      censRegLogLik <- function( beta ) {
         yMatHat <- matrix( matrix( xArr, ncol = ncol( xMat ) ) %*%
            beta[ 1:( length( beta ) - 2 ) ], nrow = nInd, ncol = nTime )
         sigmaMu <- exp( beta[ length( beta ) - 1 ] )
         sigmaNu <- exp( beta[ length( beta ) ] )
         likInd <- rep( 0, nInd )
         gradInd <- matrix( 0, nrow = nInd, ncol = length( beta ) )
         for( h in 1:nGHQ ) {
            likGhqInner <- matrix( NA, nrow = nInd, ncol = nTime )
            likGhqInner[ obsBelow ] <-
               ( left - yMatHat[ obsBelow ] - sqrt( 2 ) * sigmaMu *
                  ghqPoints$zeros[ h ] ) / sigmaNu
            likGhqInner[ obsAbove ] <-
               ( yMatHat[ obsAbove ] - right + sqrt( 2 ) * sigmaMu *
                  ghqPoints$zeros[ h ] ) / sigmaNu
            likGhqInner[ obsBetween ] <-
               ( yMat[ obsBetween ] - yMatHat[ obsBetween ] -
                  sqrt( 2 ) * sigmaMu * ghqPoints$zeros[ h ] ) / sigmaNu
            likGhq <- matrix( 1, nrow = nInd, ncol = nTime )
            likGhq[ obsBelow | obsAbove ] <-
               pnorm( likGhqInner[ obsBelow | obsAbove ] )
            likGhq[ obsBetween ] <-
               dnorm( likGhqInner[ obsBetween ] ) / sigmaNu
            likGhqProd <- apply( likGhq, 1, prod )
            likInd <- likInd + ghqPoints$weights[ h ] * likGhqProd
            # gradients
            gradPartGhq <- matrix( 0, nrow = nInd, ncol = nTime )
            gradPartGhq[ obsBelow ] <-
               - dnorm( likGhqInner[ obsBelow ] ) / sigmaNu
            gradPartGhq[ obsAbove ] <-
               dnorm( likGhqInner[ obsAbove ] ) / sigmaNu
            gradPartGhq[ obsBetween ] <-
               - ddnorm( likGhqInner[ obsBetween ] ) / sigmaNu^2
            # part of gradients with respect to beta
            for( i in 1:( length( beta ) - 2 ) ) {
               gradInd[ , i ] <- gradInd[ , i ] + ghqPoints$weights[ h ] *
                  likGhqProd * rowSums( gradPartGhq * xArr[ , , i ] / likGhq,
                  na.rm = TRUE )
            }
            # part of gradient with respect to log( sigma_mu )
            gradInd[ , length( beta ) - 1 ] <- gradInd[ , length( beta ) - 1 ] +
               sigmaMu * ghqPoints$weights[ h ] * likGhqProd *
               rowSums( gradPartGhq * sqrt( 2 ) * ghqPoints$zeros[ h ] / likGhq )
            # part of gradient with respect to log( sigma_nu )
            gradPartGhq[ obsBelow ] <- gradPartGhq[ obsBelow ] *
               likGhqInner[ obsBelow ]
            gradPartGhq[ obsAbove ] <- - gradPartGhq[ obsAbove ] *
               likGhqInner[ obsAbove ]
            gradPartGhq[ obsBetween ] <- gradPartGhq[ obsBetween ] *
               likGhqInner[ obsBetween ] - likGhq[ obsBetween ] / sigmaNu
            gradInd[ , length( beta ) ] <- gradInd[ , length( beta ) ] +
               sigmaNu * ghqPoints$weights[ h ] * likGhqProd *
               rowSums( gradPartGhq / likGhq )
         }
         ll <- log( likInd / sqrt( pi ) )
         attr( ll, "gradient" ) <- gradInd / likInd
         return( ll )
      }
   } else {
      ## naming coefficients
      names( start ) <- c( colnames( xMat ), "logSigma" )

      ## classify observations
      obsBelow <- yVec <= left
      obsAbove <- yVec >= right
      obsBetween <- !obsBelow & !obsAbove

      ## log likelihood function for cross-sectional data
      censRegLogLik <- function( beta ) {
         yHat <- xMat %*% beta[ - length( beta ) ]
         sigma <- exp( beta[ length( beta ) ] )
         ll <- rep( NA, length( yVec ) )
         ll[ obsBelow ] <-
            pnorm( ( left - yHat[ obsBelow ] ) / sigma, log.p = TRUE )
         ll[ obsBetween ] <-
            dnorm( ( yVec - yHat )[ obsBetween ] / sigma, log = TRUE ) -
            log( sigma )
         ll[ obsAbove ] <-
            pnorm( ( yHat[ obsAbove ] - right ) / sigma, log.p = TRUE )

         ## gradients of log likelihood function for cross-sectional data
         grad <- matrix( NA, nrow = length( yVec ), ncol = length( beta ) )
         grad[ obsBelow, ] <-
            dnorm( ( left - yHat[ obsBelow ] ) / sigma ) /
            pnorm( ( left - yHat[ obsBelow ] ) / sigma ) *
            cbind( - xMat[ obsBelow, , drop = FALSE ] / sigma,
               - ( left - yHat[ obsBelow ] ) / sigma )
         grad[ obsBetween, ] <-
            cbind( ( ( yVec - yHat )[ obsBetween ] / sigma ) *
               xMat[ obsBetween, , drop = FALSE ] / sigma,
               ( ( yVec - yHat )[ obsBetween ] / sigma )^2 - 1 )
         grad[ obsAbove, ] <-
            dnorm( ( yHat[ obsAbove ] - right ) / sigma ) /
            pnorm( ( yHat[ obsAbove ] - right ) / sigma ) *
            cbind( xMat[ obsAbove, , drop = FALSE ] / sigma,
               - ( yHat[ obsAbove ] - right ) / sigma )
         attr( ll, "gradient" ) <- grad
         return( ll )
      }
   }
   result <- maxLik( censRegLogLik, start = start, ... )

   # save and return the call
   result$call <- match.call()

   # save and return the number of oservations (in each category)
   result$nObs <- c( sum( obsBelow ), sum( obsBetween ), sum( obsAbove ) )
   result$nObs <- c( sum( result$nObs ), result$nObs )
   names( result$nObs ) <- c( "Total", "Left-censored", "Uncensored",
      "Right-censored" )

   class( result ) <- c( "censReg", class( result ) )
   return( result )
}

