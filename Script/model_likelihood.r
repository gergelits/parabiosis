# model_likelihood.r

calculate_model_loglik <- function( m, chain = NA )
{
  chain_dir <- get_path( "chain_dir", m = m, 
                         sel_type = "each", sel_chains = chain )
  Q <- read_Q_matrix( chain_dir = chain_dir )
  af <- read_af_matrix( m = m, chain = chain )
  model_vars_dir <- get_path( "model_vars_dir", m = m )
  load( file.path( model_vars_dir, "pb.RData" ) )
  load( file.path( model_vars_dir, "model_vars_.RData" ) )
    
  alpha <- list()
  for ( j in 1 : wn_ )
  {
    a = a0i_ %*% expm( 7 * ( w_[ j ] - w0i_ ) * Q )    
    alpha[[ j ]] <- list()
    for ( k in 1 : tn_ ) 
    {
      alpha[[ j ]][[ k ]] = c( 
        a[ ( ( k - 1 ) * pn_ + 1 ) : ( k * pn_ ) ],
        a[ ( ( tn_ + k - 1 ) * pn_ + 1 ) : ( ( tn_ + k ) * pn_ ) ] )
      alpha[[ j ]][[ k ]] = alpha[[ j ]][[ k ]] * 
        ( af[ j, k ] / sum( alpha[[ j ]][[ k ]] ) )
    }
  }

  d_ll_n <- tibble( wk = wi_, tis = ti_, ll = 0, id = 1 : n_ )
  
  loglik_parts <- tidyr::expand_grid( wk = wi_ %>% unique,
                                      tis = ti_ %>% unique,
                                      loglik = 0 ) %>% 
    arrange( wk, tis )
  
  
  # Calc log( lik_i ) for each of n_ observations, 
  # and sum separately for 3+ model parts x time points
  for ( i in 1 : n_ )
  {
    lik_i <- ( ddirichlet( x_[ i, ], alpha[[ wi_[ i ] ]][[ ti_[ i ] ]] ) )
    d_ll_n$ll[ i ] <- log( lik_i )
    
    loglik_parts$loglik[ loglik_parts$wk == wi_[ i ] & 
                           loglik_parts$tis == ti_[ i ] ] <- 
      loglik_parts$loglik[ loglik_parts$wk == wi_[ i ] & 
                             loglik_parts$tis == ti_[ i ] ] + log( lik_i )
  }
   
  # excluding Blood x wk = 0
  loglik_model_parts <- loglik_parts %>% 
    filter( ! ( wk == 1 & tis == 1 ) ) %>% 
    dplyr::rename( model_part = tis ) %>% 
    group_by( model_part ) %>% 
    summarise( ., loglik = sum( loglik ) )
  loglik <- loglik_parts %>% 
    filter( ! ( wk == 1 & tis == 1 ) ) %>% 
    summarise( sum( loglik ) ) %>% unlist()
  
  # including all parts
  loglik_model_parts_all <- loglik_parts %>% 
    dplyr::rename( model_part = tis ) %>% 
    group_by( model_part ) %>% 
    summarise( ., loglik = sum( loglik ) )
  loglik_all <- loglik_parts %>% 
    summarise( sum( loglik ) ) %>% unlist()
  
  d_ll_n_data <- cbind( d_ll_n, x_ %>% as_tibble() )
  d_ll_n_data %>% 
    arrange( tis, wk, host.celltype.naive, id ) %>% 
    write_csv( file.path( chain_dir, "ll_n_data.csv" ) )
  
  res <- list( loglik, loglik_model_parts,
               loglik_all, loglik_model_parts_all )
  names( res ) <- c( "loglik", "loglik_model_parts",
                     "loglik_all", "loglik_model_parts_all" )
  return( res )
}
