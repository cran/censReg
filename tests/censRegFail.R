library( "censReg" )

data( "Affairs", package = "AER" )

# no censored observations
try( censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs[Affairs$affairs > 10,] ) )

# no uncensored observations
try( censReg( affairs ~ age + yearsmarried + religiousness +
      occupation + rating, data = Affairs[Affairs$affairs == 0,] ) )

# standard estimations (just to get estimates)
est <- censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs )
# log-likelihood value in case of no censored observations
try( round( c( censReg( affairs ~ age + yearsmarried + religiousness +
      occupation + rating, data = Affairs[ Affairs$affairs > 10, ],
   start = coef( est ), logLikOnly = TRUE ) ), 3 ) )
# log-likelihood value in case of no censored observations
try( round( c( censReg( affairs ~ age + yearsmarried + religiousness +
      occupation + rating, data = Affairs[ Affairs$affairs == 0, ],
   start = coef( est ), logLikOnly = TRUE ) ), 3 ) )
