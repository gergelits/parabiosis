# estimate_one_markov_model.r

estimate_one_markov_model <- function( 
    m = NULL,
    use_analysis_done_as_input = FALSE,
    sel_models_file_name )
{
  ps <- m$ps
  model_name <- m$model_name
  
  # fit Markov chain model
  
  # NOTE: possible extension with ps$ANALYSIS_DONE_PATH
  model_dir <- get_path( type = "model_dir", m = m ) 
  if ( !dir.exists( model_dir ) ) { dir.create( model_dir, recursive = TRUE ) }
  
  # LOAD_PREMODEL_VARS_AND_DATA_FROM_markov_premodel_prints
  model_vars_dir <- get_path( type = "model_vars_dir", m = m )  
  load( file = sprintf( "%s/model_vars_.RData", model_vars_dir ) )
  load( file = sprintf( "%s/pb.RData", model_vars_dir ) )
  load( file = sprintf( "%s/other_vars.RData", model_vars_dir ) )
  
  # decide whether to run new MCMC simulations  
  MCMC_is_done <- file.exists( sprintf( 
    "%s/parabio_fit%i.rda", model_dir, m$mcmc_pars$mcmc.iter.n ) )                         
  if ( MCMC_is_done ) {
    load( sprintf( "%s/parabio_fit%i.rda", model_dir, m$mcmc_pars$mcmc.iter.n ) )
  } else {
    # run new MCMC simulations
    
    # copy model first to avoid access from multiple sites ( e.g., on cluster )
    file.copy( from = sprintf( "%s/%s.stan", ps$CODE_PATH, model_name ),
               to = sprintf( "%s/%s.stan", model_dir, model_name ) )
    
    model_rds <- sprintf( 
      "%s/%s/%s.rds", ps$PROJECT_PATH, "Script_model_dso", model_name )
    if ( file.exists( model_rds ) ) { 
      cat( "\nLoading exisiting, C++-compiled code ...\n\n" )
      parabio.model <- readRDS( model_rds )
    } else {
    
      cat( "\nTranslating from Stan to C++ and compiling in C++ ...\n\n" )
      
      # translate from Stan to C++ and compile in C++
      parabio.model <- stan_model( 
        file = sprintf( "%s/%s.stan", model_dir, model_name ),
        model_name = "TEST",
        warn_pedantic = TRUE, 
        warn_uninitialized = TRUE,
        save_dso = TRUE,
        verbose = FALSE
        )        
      
      # avoid rewriting at the same time
      if ( ! file.exists( model_rds ) ) {
      file.copy( from = sprintf( "%s/%s.rds", model_dir, model_name ),
                 to = model_rds ) }
    }
    
    # estimate the model by new MCMC simulation:
    cat( sprintf( "# New MCMC simulation has started: %s, %s.\n", 
                  m$celltype, m$tissue ) )
    cat( "####################################################\n\n" )
    data_list <- list( host_donor_rate_ = host_donor_rate_,
                       use_hdr2_ = use_hdr2_,
                       max_flow_to_N3_ = max_flow_to_N3_,
                       b_rate_rel_to_blood_naive_ = b_rate_rel_to_blood_naive_,                   
                       n_ = n_, wn_ = wn_, pn_ = pn_, tn_ = tn_,
                       w0i_ = w0i_, w_ = w_, wi_ = wi_, ti_ = ti_, x_ = x_,    
                       al_ = al_, a0i_ = a0i_ )
    model_sim <- purrr::quietly( sampling )(                                   
      diagnostic_file = sprintf( "%s/diagnostic_file.csv", model_dir ),
      sample_file = sprintf( "%s/sample_file.csv", model_dir ),
      verbose = TRUE,
      object = parabio.model, data = data_list, 
      iter = m$mcmc_pars$mcmc.iter.n, 
      warmup = m$mcmc_pars$mcmc.warmup.n, 
      chains = m$mcmc_pars$mcmc.chain.n,        
      control = list( adapt_delta = m$mcmc_pars$sampling.adapt_delta,
                      max_treedepth = 10 ) )                                    
    # MCMC simulation finished
    cat( sprintf( "# The MCMC simulation has finished: %s, %s.\n", 
                  m$celltype, m$tissue ) )
    cat( "####################################################\n\n" )
    
    parabio.fit <- model_sim$result
    save( parabio.fit, file = sprintf(
      "%s/parabio_fit%i.rda", model_dir, m$mcmc_pars$mcmc.iter.n ) )
    # MCMC simulation saved
    
    # save warnings, output, and messages
    sink( sprintf( "%s/parabio_fit_sampling_and_warnings.txt", model_dir ) )
    cat( "Warning messages:\n", model_sim$warnings, "\n" )
    cat( "Output:\n",           model_sim$output,   "\n" )
    cat( "Messages:\n",         model_sim$messages, "\n" )
    sink()
  }
  if ( PRINT_PARABIO_FIT_ALL_CHAINS <- TRUE ) {
    sink( file = sprintf( "%s/parabio_fit_print.txt", model_dir ) )
    width.backup <- getOption( "width" ); options( width = 200 ) 
    print( parabio.fit )
    options( width = width.backup )
    sink()
  }
  
  
  # Which chains to get results for:

  sel_chains_table <- tibble( 
    sel_chains_set = c( 1 : m$mcmc_pars$mcmc.chain.n ) %>% as.character(),
    sel_chains_type = c( rep( "each", m$mcmc_pars$mcmc.chain.n ) ) )
  write_csv( sel_chains_table, sprintf( "%s/sel_chains_table_TMP.csv", model_dir ) )
  
  # either all chain, which takes ~5 minutes:
  for( i in 1 : nrow( sel_chains_table ) ) {
    
  # or only max lp__ chain to save time:
  # for( i in ( which( sel_chains_table$sel_chains_type == "max" ) ) )   
    
    get_chain_output( 
      m = m, parabio.fit = parabio.fit,
      sel_chains_table = sel_chains_table, i = i, 
      sel_models_file_name = sel_models_file_name )
  }

  # Add 3mad chains calculation  
  sel_chains_set_3mad <- select_mcmc_chains(
    parabio.fit = parabio.fit, 
    m = m,
    mcmc.chain.n = m$mcmc_pars$mcmc.chain.n, 
    type = "3mad" ) %>%
    paste( collapse = "_" )
  if ( sel_chains_set_3mad != "" ) {
    sel_chains_table <- bind_rows( sel_chains_table, 
                                   tibble( sel_chains_set = sel_chains_set_3mad, 
                                           sel_chains_type = "3mad" ) )
  }
  get_chain_output( 
    m = m, parabio.fit = parabio.fit,
    sel_chains_table = sel_chains_table, i = nrow( sel_chains_table ), 
    sel_models_file_name = sel_models_file_name )
  
  # Add max chain calculation  
  sel_chains_set_max <- select_mcmc_chains( 
    parabio.fit = parabio.fit, 
    m = m,
    mcmc.chain.n = m$mcmc_pars$mcmc.chain.n, 
    type = "max" ) %>%
    paste( collapse = "_" )
  sel_chains_table <- bind_rows( sel_chains_table, 
                                 tibble( sel_chains_set = sel_chains_set_max,
                                         sel_chains_type = "max" ) )
  get_chain_output( 
    m = m, parabio.fit = parabio.fit,
    sel_chains_table = sel_chains_table, i = nrow( sel_chains_table ), 
    sel_models_file_name = sel_models_file_name )
  
  write_csv( sel_chains_table, sprintf( "%s/sel_chains_table.csv", model_dir ) )
}
