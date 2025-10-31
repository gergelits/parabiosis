# calculate_parabiosis_paper.r
 
# Copyright (c) 2025 University of Cambridge and Babraham Institute (United Kingdom)
#
# Software written by Vaclav Gergelits and Carlos P. Roca as research funded 
# by the European Union and the United Kingdom
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.

# Runs MCMC simulation, simulation of cell population based on Markov chain model,
# and creates all figures and tables used in the manuscript.


calculate_parabiosis_paper <- function( 
    ps = PS,   
    CELLTYPE = "Treg",
    # the tissue modelled by Markov Chain
    F_TISSUES = NA,
    # other tissues standing in the model separately
    VECTOR_TISSUES_0 = NULL,
    # use precalculated MCMC simulations
    use_analysis_done_as_input = FALSE,
    # stan model file name (without extension .stan)
    MODEL_NAME = "model_m0003",
    # stan_pars
    stan_pars_v = "v01",
    # mcmc_pars
    mcmc_pars_v = "i1001",
    # combined dwell time simulations
    n_sim = 1000 )

{
  # A few Tissue - Cell type combinations that are not calculated
  if ( F_TISSUES == "Thymus" ) { return( NULL ) }
  if ( CELLTYPE %in% c( "NK", "B" ) & ( F_TISSUES == "BoneMarrow" ) ) {
    return( NULL )
  }
  
  # load libraries and my own functions
  source( sprintf( "%s/%s", ps$CODE_PATH, "install_packages.r" ) )
  source( sprintf( "%s/%s", ps$CODE_PATH, "load_packages.r" ) )
  load_packages()
  source( sprintf( "%s/%s", ps$CODE_PATH, "load_my_functions.r" ) )
  load_my_functions( ps = ps )
  
  
  # create directories, get parameters
  model_id <- get_model_id( file_name = MODEL_NAME, from_file = "stan" )
  if ( VECTOR_TISSUES_0 == "NULL" ) VECTOR_TISSUES_0 <- NULL
  create_dirs( ps = ps )
  st <- get_stan_pars( stan_pars_v = stan_pars_v, ps = ps )
  mcmc_pars <- get_mcmc_pars( mcmc_pars_v = mcmc_pars_v, ps = ps )  
  
  
  # setup modelling parameters
  if ( ! is.na( F_TISSUES ) ) {
    # tissue argument passed -> compute only 1 tissue
    TISSUES <- F_TISSUES 
  } else { 
    # tissue argument not passed -> compute all tissues
    TISSUES <- get_TISSUES( celltype = CELLTYPE )
  } 
  N_TISSUES <- length( VECTOR_TISSUES_0 ) + 1
  st$MODEL_VER <- sprintf( "%s_t%s", stan_pars_v, N_TISSUES )
  # NOTE:
  # Possible extension: extra comment for a small change - in anything
  # if ( !is.na( st$F_MODEL_VER_EXT ) ) { 
  #   st$MODEL_VER = sprintf( "%s%s", st$MODEL_VER, st$F_MODEL_VER_EXT ) }      
  # If extended then write get_MODEL_VER function
  m_tisNA <- model( ps = ps, celltype = CELLTYPE,
                    mcmc_pars = mcmc_pars, stan_pars_v = stan_pars_v, st = st,
                    model_name = MODEL_NAME, model_ver = st$MODEL_VER )
  SEL_MODELS_FILE_NAME <- get_sel_models_file_name( m_tisNA )

  
  # for each defined tissue calculate MCMC model 
  # this is only one tissue when run from Shell script
  for ( TIS_I in TISSUES )
  {
    m <- model( m_tisNA, tissue = TIS_I )
    VECTOR_TISSUES <- c( TIS_I, VECTOR_TISSUES_0 )
    if( length( VECTOR_TISSUES ) > length( unique( VECTOR_TISSUES ) ) )
      { VECTOR_TISSUES <- c( unique( VECTOR_TISSUES ), "LN" ) }
    COMPLEM_TISSUES = setdiff( 
      get_TISSUES( celltype = CELLTYPE ), VECTOR_TISSUES )
    # Compute new MCMC simulation  
    # TISSUES ............ typically only 1 tissue (backward compatibility)
    # TIS_I .............. that only tissue from TISSUES
    # VECTOR_TISSUES ..... all tissues standing separately in the MC model
    # VECTOR_TISSUES_0 ... other standing separately in the MC model
    # COMPLEM_TISSUES .... all other tissues in the model - are pooled

    markov_premodel_prints_eps( 
      m = m,
      vector_tissues = VECTOR_TISSUES,
      complem_tissues = COMPLEM_TISSUES ) 
    estimate_one_markov_model(
      m = m,
      use_analysis_done_as_input = use_analysis_done_as_input,
      sel_models_file_name = SEL_MODELS_FILE_NAME ) 
    while ( sink.number() > 0 ) sink()
    
    # calculate the MCMC accuracy:  
    m <- model( m = m, max_lp_chain = get_max_chain( m ) )
    calculate_HDI(
      m = m, chain = "max",
      use_analysis_done_as_input = use_analysis_done_as_input )
    
    write_model_lp_loglik( m = m )
  }
  
  
  # plot figures for all the calculated tissues
  write_sel_models_file(
    m = m, celltype = CELLTYPE, ps = ps, mcmc_pars = mcmc_pars, 
    model_name = MODEL_NAME, model_ver = st$MODEL_VER,
    sel_models_file_name = SEL_MODELS_FILE_NAME )
  m_tisNA <- get_m_tisNA( m = m )
  visualize_results( m_tisNA = m_tisNA )  
  
  
  for ( TIS_I in TISSUES ) {
    # Figure only for a specific tissue at the end here - it takes extra time
    sim_10000 <- sim_combined_dwell( m = m, n_sim = n_sim )
    plot_combined_dwell( m = m, sim = sim_10000, export = TRUE )
    plot_16_comb_dwells( m_tisNA = get_m_tisNA( m = m ), n_sim = n_sim )
  }  
  plot_exit_types( m_tisNA = m_tisNA, n_sim = n_sim )
  get_rel_flows_and_exits( m_tisNA = m_tisNA, n_sim = n_sim )
  plot_influx( m = m_tisNA, n_sim = n_sim )
  
  plot_all_celltypes_dwells( m_tisNA = m_tisNA, n_sim = n_sim )
  plot_fig_7_comb_dwells( m_tisNA = m_tisNA, n_sim = n_sim, ylim_max = 50 )
  
  if ( TRUE ) { sessionInfo() }
  cat( "calculate_parabiosis_paper() finished.\n\n" )
}
