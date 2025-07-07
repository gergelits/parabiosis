# setup_simulation_parameters_i101.r

seed.base <- 20190601
seed <- set_seed_here_visible( seed.base, "markov model" )
mcmc.iter.n <- 200
mcmc.warmup.n <- 100
mcmc.chain.n <- 1
mcmc.core.n <- mcmc.chain.n
sampling.adapt_delta <- 0.50
options( mc.cores = mcmc.core.n )
rstan_options( auto_write = TRUE )
