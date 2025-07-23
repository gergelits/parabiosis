# get_mcmc_pars.r

get_mcmc_pars <- function( mcmc_pars_v, ps = ps )
{
  mcmc_pars_script <- sprintf( "setup_simulation_parameters_%s.r", mcmc_pars_v )
  source( sprintf( "%s/%s", ps$CODE_PATH, mcmc_pars_script ), local = TRUE )
  
  mcmc_pars <- list()
  mcmc_pars$seed <- seed
  mcmc_pars$mcmc.iter.n <- mcmc.iter.n
  mcmc_pars$mcmc.warmup.n <- mcmc.warmup.n
  mcmc_pars$mcmc.chain.n <- mcmc.chain.n
  mcmc_pars$sampling.adapt_delta <- sampling.adapt_delta
  mcmc_pars$n.iter <- mcmc.chain.n * ( mcmc.iter.n - mcmc.warmup.n )
  mcmc_pars$id <- mcmc_pars_v
  # mcmc_pars_v must be coding number of iterations after warmup and n of chains, 
  # e.g., 5004 -> 5000 iterations after warmup on each of 4 chains
  mcmc_pars_v_i <- substring( mcmc_pars$id, first = 2 ) %>% as.numeric
  stopifnot( mcmc_pars_v_i %/% 100 == ( mcmc.iter.n - mcmc.warmup.n ) %/% 100 )
  stopifnot( mcmc_pars_v_i %% 100 == mcmc_pars$mcmc.chain.n )
  
  return( mcmc_pars )
}
