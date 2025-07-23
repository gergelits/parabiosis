# setup_simulation_parameters_i5004.r

seed.base <- 20230101
seed <- set_seed_here_visible( seed.base, "markov model" )
mcmc.iter.n <- 5500
mcmc.warmup.n <- 500
mcmc.chain.n <- 4
mcmc.core.n <- mcmc.chain.n
sampling.adapt_delta <- 0.50
options( mc.cores = mcmc.core.n )
rstan_options( auto_write = TRUE )
