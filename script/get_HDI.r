#  get_HDI.r

get_HDI_4_values <- function( MCMC.sample, 
                              credMass.big = 0.80, credMass.small = 0.5 )
{
  credMass <- NULL
  sorted_sample <- sort( MCMC.sample )
  ciIdxInc_b <- ceiling( credMass.big * length( sorted_sample ) )
  ciIdxInc_s <- ceiling( credMass.small * length( sorted_sample ) )
  nCIs_b <- length( sorted_sample ) - ciIdxInc_b
  nCIs_s <- length( sorted_sample ) - ciIdxInc_s
  ciWidth_b <- rep( 0, nCIs_b )
  ciWidth_s <- rep( 0, nCIs_s )
  for ( i in 1 : nCIs_b ) {
    ciWidth_b[ i ] <- sorted_sample[ i + ciIdxInc_b ] - sorted_sample[ i ]
  }
  for ( i in 1 : nCIs_s ) {
    ciWidth_s[ i ] <- sorted_sample[ i + ciIdxInc_s ] - sorted_sample[ i ]
  }
  HDImin_b <- sorted_sample[ which.min( ciWidth_b ) ]
  HDImax_b <- sorted_sample[ which.min( ciWidth_b ) + ciIdxInc_b ]
  HDImin_s <- sorted_sample[ which.min( ciWidth_s ) ]
  HDImax_s <- sorted_sample[ which.min( ciWidth_s ) + ciIdxInc_s ]
  HDIlim_4 <- sort( c( HDImin_b, HDImin_s, HDImax_s, HDImax_b ) )
  return( HDIlim_4 )
}
