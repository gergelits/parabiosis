# select_mcmc_chains.r

select_mcmc_chains <- function ( parabio.fit, m = "TODO", mcmc.chain.n, type = "3mad" )
{
  fit.lp <- rep( NA, mcmc.chain.n )
  fit.ll <- rep( NA, mcmc.chain.n )
  for ( i in 1 : mcmc.chain.n ) {
    try( fit.lp[ i ] <- get_posterior_mean( parabio.fit )[ "lp__", i ] )
    try( fit.ll[ i ] <- 
           calculate_model_loglik( m = m, chain = i )$loglik )
  }
  
  if ( type == "3mad" )
  {
    fit.lp <- ifelse( fit.lp < 0, NA, fit.lp )
    
    if ( sum( !is.na( fit.lp ) ) == 0 )
    {
      parabio.fit.chain <- c()
    } 
    
    if ( sum( !is.na( fit.lp ) ) == 1 )
    {
      parabio.fit.chain <- which( !is.na( fit.lp ) )
    } 
    
    if ( sum( !is.na( fit.lp ) ) >= 2 )
    {
      parabio.fit.chain <- ( 1 : mcmc.chain.n )[
        which( abs( fit.lp - median( fit.lp, na.rm = TRUE) ) / 
                 stats::mad( fit.lp, na.rm = TRUE ) < 3 ) ]
    }
  }
  
  if ( type == "max" ) {
    parabio.fit.chain <- which.max( fit.ll )
  }
  
  return( parabio.fit.chain )
}


check_MCMC_convergence <- function( m = NULL ) {
  model_dir <- get_path( "model_dir", m = m )
  MCMC_convergence <- file.exists( 
    sprintf( "%s/parabio_fit%i.rda", model_dir, m$mcmc_pars$mcmc.iter.n ) ) & 
    ( file.exists( sprintf( "%s/parabio_fit_print.txt", model_dir ) ) ) & 
    ( file.info( sprintf( "%s/parabio_fit_print.txt", model_dir ) )$size > 2000 )
  return( MCMC_convergence )
}


write_model_lp_loglik <- function( m )
{
  # NOTE: 
  # rstan::get_logposterior()
  # Get the log-posterior at each iteration. Each element of the returned list 
  # is the vector of log-posterior values (up to an additive constant, i.e. up 
  # to a multiplicative constant on the linear scale) for a single chain. 
  # The optional argument inc_warmup (defaulting to TRUE) indicates whether to 
  # include the warmup period.
  # rstan::log_prob()
  # Compute the log probability density (lp__) for a set of parameter values 
  # (on the unconstrained space) up to an additive constant. The unconstrained 
  # parameters are specified using a numeric vector. The number of parameters 
  # on the unconstrained space can be obtained using method get_num_upars. 
  # A numeric value is returned. See also the documentation in log_prob.
  
  model_dir <- get_path( "model_dir", m = m )
  sel_chains_table <- read_csv( file.path( model_dir,
                                           "sel_chains_table.csv" ) ) %>% 
    # compatibility with older versions:
    dplyr::select( sel_chains_set, sel_chains_type ) %>% 
    dplyr::rename( chain = sel_chains_set )
  
  ch_diag <- 
    sel_chains_table %>% 
    # NOTE: Possible to extend. sel_chains_type == ...
    filter( sel_chains_type == "each" ) %>%         
    mutate( lp__mean = NA,
            lp__mode = NA,
            loglik = NA, 
            loglik_all = NA,
            loglik_p1 = NA, loglik_p1_all = NA, loglik_p2 = NA, loglik_p3 = NA,
            t_sample_s = NA, t_warmup_s = NA )
  stopifnot( m$mcmc_pars$mcmc.chain.n == ( ch_diag$chain %>% max ) )
  
  for ( i in 1 : nrow( ch_diag ) )
  {
    chain <- ch_diag$chain[ i ]
    parabio_fit_HDI <- get_parabio_fit_HDI( m = m, chain = chain )
    ch_diag$lp__mean[ i ] <- 
      parabio_fit_HDI$mean[ parabio_fit_HDI$par == "lp__" ]
    ch_diag$lp__mode[ i ] <- 
      parabio_fit_HDI$mode[ parabio_fit_HDI$par == "lp__" ]
    loglik <- calculate_model_loglik( m = m, chain = chain )
    ch_diag$loglik[ i ] <- loglik$loglik
    ch_diag$loglik_all[ i ] <- loglik$loglik_all
    ch_diag$loglik_p1[ i ] <- loglik$loglik_model_parts$loglik[ 1 ]
    ch_diag$loglik_p1_all[ i ] <- loglik$loglik_model_parts_all$loglik[ 1 ]
    ch_diag$loglik_p2[ i ] <- loglik$loglik_model_parts$loglik[ 2 ]
    ch_diag$loglik_p3[ i ] <- loglik$loglik_model_parts$loglik[ 3 ]
  }
  
  chains_elapsed_time_ff <- file.path( model_dir, "chains_elapsed_time.csv" )
  if ( file.exists( chains_elapsed_time_ff ) ) {
    table_elapsed_time <- 
      read_csv( chains_elapsed_time_ff ) %>% 
      mutate( chain = as.character( chain ) )
    ch_diag <- 
      ch_diag %>% 
      dplyr::select( -c( t_sample_s, t_warmup_s ) ) %>%
      mutate( chain = as.character( chain ) ) %>% 
      left_join( x = ., y = table_elapsed_time, by = c( "chain" = "chain" ) )
  }
  write_csv( ch_diag, file.path( model_dir, "each_chain_diag_table.csv" ) )
  
  figs_diag_dir <- get_path( "figs_m_diag", m = m )
  if ( !dir.exists( figs_diag_dir ) ) { dir.create( figs_diag_dir, 
                                                    recursive = TRUE ) }
  chain_dir <- get_path( 
    "chain_dir", m = m, sel_type = "max", sel_chains = get_max_chain( m ) )
  
  file.copy(
    from = file.path( model_dir, "model_vars", "al_.csv" ),
    to = file.path( figs_diag_dir, sprintf( "al__%s.csv", m$tissue ) ) )
  file.copy(
    from = file.path( chain_dir, "Q_matrix.csv" ),
    to = file.path( figs_diag_dir, sprintf( "Q_%s_maxlp.csv", m$tissue ) ) )
  file.copy(
    from = file.path( model_dir, "each_chain_diag_table.csv" ),
    to = file.path( figs_diag_dir, sprintf( 
      "each_chain_diag_table_%s.csv", m$tissue ) ) )
  
  chains_elapsed_time_ff
  if ( file.exists( chains_elapsed_time_ff ) ) { ch_diag }
}
