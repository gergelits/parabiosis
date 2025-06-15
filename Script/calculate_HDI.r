# calculate_HDI.r


calculate_HDI <- function( m = NULL, chain = "max",
                           use_analysis_done_as_input = FALSE )     
{
  # NOTE: possible extension with ps$ANALYSIS_DONE_PATH
  to_calc_HDI <- TRUE
  if ( ! check_MCMC_convergence( m = m ) ) return( NULL ) 
  
  # Now, successful MCMC simulation exists
  model_dir <- get_path( "model_dir", m = m )
  load( sprintf( "%s/parabio_fit%i.rda", model_dir, m$mcmc_pars$mcmc.iter.n ) )
  if ( ! file.exists( file.path( model_dir, "chains_elapsed_time.csv" ) ) )
  {
    try( {
      elapsed_time <- rstan::get_elapsed_time( parabio.fit )
      elapsed_time %>% 
        as_tibble %>% 
        mutate( chain = 1 : nrow( . ) ) %>% 
        dplyr::rename( t_warmup_s = warmup,
                       t_sample_s = sample ) %>% 
        dplyr::select( chain, everything() ) %>% 
        write_csv( file.path( model_dir, "chains_elapsed_time.csv" ) )
    } )
  }
    
  if ( chain == "max" )
  {
    parabio.fit.chain <- get_max_chain( m )
    out_dir <- model_dir
  } else {                           
    # NOTE: Possible extension: check 3mad
    # chain is an integer
    stopifnot( chain %in% 1 : 20 )
    parabio.fit.chain <- chain
    out_dir <- get_path( "chain_dir", m = m, 
                         sel_type = "each", sel_chains = chain )
  }
  if ( ! dir.exists( out_dir ) ) { dir.create( out_dir ) }
  
  
  if ( file.exists( sprintf( "%s/parabio_fit_HDI.csv", out_dir ) ) )
  {
    to_calc_HDI <- FALSE
    cat( sprintf( "\nparabio_fit_HDI.csv already exists for %s, chain %s.\n", 
                  m$tissue, chain ) )
  } 
    
  to_calc_HDI <- TRUE
  # from this point general for any parabio.fit.chain 
  if ( to_calc_HDI ) {
    parabio.fit.sample <- rstan::extract( parabio.fit, permuted = FALSE )
    dParamsHDI <- NULL
    for ( mp in dimnames( parabio.fit.sample )$parameters )
    {
      parabio.fit.chain <- as.numeric( parabio.fit.chain )
      pf.sample <- as.vector( parabio.fit.sample[ , parabio.fit.chain, mp ] )
      pf.sample.mode <- mlv( pf.sample, method = "venter", type = "shorth" )
      HDI.80.50 <- get_HDI_4_values( MCMC.sample = pf.sample, 
                                     credMass.big = 0.80, 
                                     credMass.small = 0.5 )
      dParamsHDI <- 
        dParamsHDI %>% 
        bind_rows( tibble( 
          par = mp, 
          mode = pf.sample.mode, 
          hdi.80.lower = HDI.80.50[ 1 ] %>% as.numeric, 
          hdi.80.upper = HDI.80.50[ 4 ] %>% as.numeric,
          hdi.50.lower = HDI.80.50[ 2 ] %>% as.numeric, 
          hdi.50.upper = HDI.80.50[ 3 ] %>% as.numeric,
          mean = mean( pf.sample, na.rm = TRUE ),
          q50 =   quantile( pf.sample, 0.500 ),
          q2.5 =  quantile( pf.sample, 0.025 ),
          q97.5 = quantile( pf.sample, 0.975 ),
          q10 =   quantile( pf.sample, 0.100 ),
          q90 =   quantile( pf.sample, 0.900 ),
          q25 =   quantile( pf.sample, 0.250 ),
          q75 =   quantile( pf.sample, 0.750 ) ) )
    }
    dParamsHDI <- 
      dParamsHDI %>%
      mutate( 
        hdi.80.ratio = hdi.80.upper / hdi.80.lower,
        hdi.50.ratio = hdi.50.upper / hdi.50.lower )
    
    dParamsHDI %>% 
      write_csv( ., sprintf( "%s/parabio_fit_HDI.csv", out_dir ) ) 
  }
  return( TRUE )
}
