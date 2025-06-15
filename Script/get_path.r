# get_path.r

not_null_na <- function( obj ) { ! is.null( obj ) && ! is.na( obj ) }

model <- function( model = NULL, m = NULL,
                   ps = NA, celltype = NA, tissue = NA, 
                   model_name = NA, model_ver = NA, 
                   mcmc_pars = NA, 
                   stan_pars_v = NA, st = NA,
                   max_lp_chain = NA, chain = NA )
{
  if ( is.null( m ) ) { m <- model }
  
  if( is.null( m ) ) { 
    m <- list() 
    m$ps <- ps
    m$celltype <- celltype
    m$tissue <- tissue
    m$model_name <- model_name
    m$model_ver <- model_ver
    m$mcmc_pars <- mcmc_pars
    m$stan_pars_v <- stan_pars_v
    m$st <- st
    m$max_lp_chain <- max_lp_chain
    m$chain <- chain
  } else {
    if ( not_null_na( ps ) ) { m$ps <- ps }
    if ( ! is.na( celltype ) ) { m$celltype <- celltype }
    if ( ! is.na( tissue ) ) { m$tissue <- tissue }
    if ( ! is.na( model_name ) ) { m$model_name <- model_name }
    if ( ! is.na( model_ver ) ) { m$model_ver <- model_ver }
    if ( not_null_na( mcmc_pars ) ) { m$mcmc_pars <- mcmc_pars }
    if ( ! is.na( stan_pars_v ) ) { m$stan_pars_v <- stan_pars_v }
    if ( not_null_na( st ) ) { m$st <- st }
    if ( not_null_na( max_lp_chain ) ) { m$max_lp_chain <- max_lp_chain }
    if ( not_null_na( chain ) ) { m$chain <- chain }
  }
  if ( ! is.na( m$model_name ) & is.null( m$model_id ) ) { 
    m$model_id <- substring( 
      model_name, 
      first = nchar( "model_" ) + 1, last = nchar( "model_" ) + 5 ) 
  }
  return( m )
}


get_m_tisNA <- function( m )
{
  m_tisNA <- m
  m_tisNA$tissue <- NULL
  m_tisNA$max_lp_chain <- NULL
  return( m_tisNA )
}



get_fig_fname_start <- function( model = NULL, m = NULL, celltype, 
                                 model_id, mcmc_pars, combined = FALSE )
{
  if ( is.null( m ) ) { m <- model }
  if ( combined ) { m$celltype <- "Combined" }
  if ( ! is.null( m ) ) {
    start <- sprintf( "%s_%s_%s_%s", 
                      m$celltype, m$model_id, m$model_ver, m$mcmc_pars$id )
  } else {
    start <- sprintf( "%s_%s_%s", celltype, model_id, mcmc_pars$id )
  }
  return( start )
}


get_path <- function( 
    type = "MISSING" , 
    model = NULL, m = NULL,
    celltype, tissue, ps, model_dir = NA, model_dir_r = NA,
    model_name, model_id, mcmc_pars, model_ver,
    sel_type, sel_chains, chain )
{
  if ( is.null( m ) ) { m <- model }
  if( ! is.null( m ) ) {
    if ( ! is.null( m$ps ) )         { ps <- m$ps }
    if ( ! is.na( m$celltype ) )     { celltype <- m$celltype }
    if ( not_null_na( m$tissue ) )   { tissue <- m$tissue }
    if ( ! is.na( m$model_name ) )   { model_name <- m$model_name }
    if ( ! is.na( m$model_ver ) )    { model_ver <- m$model_ver }
    if ( ! is.null( m$mcmc_pars ) )  { mcmc_pars <- m$mcmc_pars }
    if ( ! is.na( m$stan_pars_v ) )  { stan_pars_v <- m$stan_pars_v }
    if ( ! is.null( m$st ) )         { st <- m$st }
    if ( not_null_na( m$max_lp_chain ) ) { max_lp_chain <- m$max_lp_chain }
    if ( not_null_na( m$model_id ) ) { model_id <- m$model_id }  
  }
  
  if ( ( ! exists( "model_id" ) ) & ( exists( "m" ) ) ) { 
    model_id <- get_model_id( m ) }
  
  if ( type == "analysis_tissue_dir" )
  {
    path <- sprintf( "%s/%s/%s", ps$ANALYSIS_PATH, celltype, tissue )
  }
  
  # r = relative path
  if ( type == "model_dir_r" ) 
  {
    if ( ! exists( "model_id" ) ) { 
      model_id <- substring( model_name, 
                             first = nchar( "model_" ) + 1, 
                             last = nchar( "model_" ) + 5 ) }
    if ( ! is.na( model_ver ) ) {
      model_dir_r <- sprintf( "%s_%s_%s", model_id, model_ver, mcmc_pars$id )
    } else {
      model_dir_r <- sprintf( "%s_%s", model_id, mcmc_pars$id )
    }
    path <- model_dir_r
  }
  
  if ( type == "model_dir" ) 
  {
    analysis_tissue_dir <- get_path( 
      "analysis_tissue_dir", m = m,
      ps = ps, celltype = celltype, tissue = tissue )
    if ( is.na( model_dir_r ) ) {
      if ( ! exists( "model_id" ) ) { 
        model_id <- substring( model_name, 
                               first = nchar( "model_" ) + 1, 
                               last = nchar( "model_" ) + 5 ) }
      
      model_dir_r <- get_path( 
        "model_dir_r", m = m,
        model_id = model_id, model_ver = m$model_ver,
        mcmc_pars = mcmc_pars )
    }
    path <- sprintf( "%s/%s", analysis_tissue_dir, model_dir_r )
  }
  
  if ( type == "model_vars_dir" )
  {
    if ( ! exists( "model_id" ) ) { 
      model_id <- substring( model_name, 
                             first = nchar( "model_" ) + 1, 
                             last = nchar( "model_" ) + 5 ) }
    path <- sprintf( 
      "%s/%s",
      get_path( "model_dir", m = m,
                ps = ps, celltype = celltype, tissue = tissue,
                model_id = model_id, model_ver = model_ver,
                mcmc_pars = mcmc_pars ),
      "model_vars" )
  }
  
  if ( type == "chain_dir" )
  {
    if ( ! exists( "model_id" ) ) { 
      model_id <- substring( model_name, 
                             first = nchar( "model_" ) + 1, 
                             last = nchar( "model_" ) + 5 ) }
    if ( is.na( model_dir ) ) {
      model_dir <- get_path( 
        "model_dir", m = m,
        ps = ps, celltype = celltype, tissue = tissue,
        model_id = model_id, model_ver = model_ver,
        mcmc_pars = mcmc_pars )
    }
    path <- sprintf( "%s/%s_chains_%s", model_dir, sel_type, sel_chains )
  }

  if ( type == "data" )
  {
    path <- sprintf( "%s/%s", ps$PROCESSED_PATH, celltype )
  }
  
  if ( type == "results" )
  {
    path <- sprintf( "%s/%s", ps$RESULTS_PATH, celltype )   
  }
  
  if ( type == "figs_comb" )
  {
    path <- sprintf( "%s/%s", ps$RESULTS_PATH, "Combined" )   
  }
  
  if ( type == "figs" )
  {
    path <- sprintf( "%s/%s/%s", ps$RESULTS_PATH, celltype, "Figures" ) 
  }
  
  if ( type == "figs_m" )
  {
    if ( ! is.null( m ) )
    {
      path <- sprintf( "%s/%s/%s/%s", ps$RESULTS_PATH, celltype, "Figures", 
                       get_fig_fname_start( m ) ) 
    } else {
      path <- sprintf( "%s/%s/%s", ps$RESULTS_PATH, celltype, "Figures" ) 
    }
  }

  if ( type == "figs_m_traj" )
  {
    if ( ! is.null( m ) )
    {
      path <- sprintf( "%s/%s/%s/%s/%s", ps$RESULTS_PATH, celltype, "Figures", 
                       get_fig_fname_start( m ), "traj" ) 
    } else {
      path <- sprintf( "%s/%s/%s/%s", 
                       ps$RESULTS_PATH, celltype, "Figures", "traj" ) 
    }
  }
    
  if ( type == "figs_m_dwell" )
  {
    if ( ! is.null( m ) )
    {
      path <- sprintf( "%s/%s/%s/%s/%s", ps$RESULTS_PATH, celltype, "Figures", 
                       get_fig_fname_start( m ), "dwell" ) 
    } else {
      path <- sprintf( "%s/%s/%s/%s", 
                       ps$RESULTS_PATH, celltype, "Figures", "dwell" ) 
    }
  }
  
  if ( type == "figs_m_rda" )
  {
    if ( ! is.null( m ) )
    {
      path <- sprintf( "%s/%s/%s/%s/%s", ps$RESULTS_PATH, celltype, "Figures", 
                       get_fig_fname_start( m ), "rda" ) 
    } else {
      path <- sprintf( "%s/%s/%s/%s", 
                       ps$RESULTS_PATH, celltype, "Figures", "rda" ) 
    }
  }
  
  if ( type == "figs_m_diag" )
  {
    if ( ! is.null( m ) )
    {
      path <- sprintf( "%s/%s/%s/%s/%s", ps$RESULTS_PATH, celltype, "Figures", 
                       get_fig_fname_start( m ), "diag" ) 
    } else {
      path <- sprintf( "%s/%s/%s/%s", 
                       ps$RESULTS_PATH, celltype, "Figures", "diag" ) 
    }
  }
  
  if ( type == "flow_fig" )
  {
    path <- sprintf( "%s/%s/Flow_Diagram", ps$RESULTS_PATH, celltype )  
  }

  return( path )
}



get_figs_rda_name <- function( fig_filename_start, 
                               end = "gg_figs_all_2.rda" )
{
  name <- sprintf( "%s_%s", fig_filename_start, end )
  return( name )
}


get_model_id <- function( file_name, from_file = "sel_models", 
                          model = NULL, m = NULL )
{
  if ( is.null( m ) ) { m <- model }
  
  if ( ! is.null( m ) & ! is.null( m$model_name ) ) {
    model_id <- substring( m$model_name, 
                           first = nchar( "model_" ) + 1, 
                           last = nchar( "model_" ) + 5 )
  } else {
    if ( from_file == "sel_models" )
    {
      model_id <- sprintf( "m%s", 
                           gsub( pattern = ".*_m([0-9]{3,4}[a-z]{0,3})_.*", 
                           replacement = "\\1", x = file_name ) )
    }
    
    if ( from_file == "stan" )
    {
      model_id <- substring( 
        file_name, first = nchar( "model_" ) + 1, last = nchar( "model_" ) + 5 )
    }
  }
  return( model_id )
}


get_max_chain <- function( m )
{
  model_dir <- get_path( type = "model_dir", m = m )
  sel_chains_csv <- sprintf( "%s/sel_chains_table.csv", model_dir )
  if ( file.exists( sel_chains_csv ) ) {
    sel_chains_table <- read.csv( sel_chains_csv )
    max_chain <- sel_chains_table$sel_chains_set[ 
      sel_chains_table$sel_chains_type == "max" ]
    
    # older version:
    if ( is.null( max_chain ) ) {
      max_chain <- sel_chains_table$sel.chains_set[ 
        sel_chains_table$sel.chains_type == "max" ]
    }
  } else {
    max_chain <- NA
  }
  return( max_chain )
}


get_sel_models_file_name <- function( m_tisNA )
{
  name <- sprintf(
    "%s_sel_models_%s_%s_%s.csv",
    m_tisNA$celltype,
    substring( m_tisNA$model_name, first = nchar( "model_" ) + 1, 
               last = nchar( "model_" ) + 5 ),
    m_tisNA$stan_pars_v,
    m_tisNA$mcmc_pars$id )
  return( name )
}


read_load_multi <- function( fun = "read_csv", full_path )
{
  f <- eval( parse( text = fun ) )
  op <- options( digits.secs = 6 )
  set_seed_here_visible( 0, Sys.time() ) 
  options( op )
  
  full_path_tmp <- sprintf( "%s_%i", full_path, sample.int( 100000, size = 1 ) )
  file.copy( from = full_path, to = full_path_tmp )
  
  # to correct when the file is not copied
  i_safety <- 0
  while ( file.size( full_path_tmp ) <= 1 & ( i_safety <= 10 ) ) {
    cat( "correction - empty file", i_safety, "\n" )
    Sys.sleep( runif( 1, min = i_safety / 2 + 0.1, max = i_safety + 0.1 ) )
    file.remove( full_path_tmp )
    full_path_tmp <- sprintf( "%s_%i", full_path, 
                              sample.int( 100000, size = 1 ) )
    file.copy( from = full_path, to = full_path_tmp )
    i_safety <- i_safety + 1
  }
  
  file <- f( full_path_tmp )
  file.remove( full_path_tmp )
   
  return( file )
}
