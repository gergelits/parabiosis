# converting_Q_and_flow_tables.r

get_exit_flow_table <- function( ps, celltype, tissue, model_dir_r )
{
  model_dir <- get_path( "model_dir", m = NULL,
                         ps = ps, celltype = celltype, tissue = tissue,
                         model_dir_r = model_dir_r )
  exit.flow.file <- sprintf( "%s/parabio_fit_HDI.csv", model_dir )
  
  if ( file.exists( exit.flow.file ) ) {
    exit.flow.table <- 
      read_csv( exit.flow.file, show_col_types = FALSE ) %>% 
      dplyr::select( par, mode ) %>% 
      filter( par %in% c( "q12", "q23", "q32" ) |
        stringr::str_detect( par, "^q[1-6][1-6]_i\\[1" ) ) %>% 
      filter( ! ( par %in% sprintf( "q%s%s_i[1]", 1:6, 1:6 ) ) ) %>%
      mutate( par = substring( par, first = 1, last = 3 ) ) %>% 
      
      mutate( population.1 = str_match( par, "^q(\\d)[1-6]" )[ , 2 ] %>% 
                as.integer,
              population.2 = str_match( par, "^q[1-6](\\d)" )[ , 2 ] %>% 
                as.integer ) %>% 
      dplyr::rename( flow = mode ) %>% 
      mutate( flow.part = "exit" ) %>% 
      dplyr::select( par, population.1, population.2, flow, flow.part ) %>% 
      arrange( population.1, population.2 )
    return( exit.flow.table )
  }
}


get_al_ <- function( model = NULL, m = NULL,
                     ps, celltype, tissue, model_dir_r = NA,
                     # if model_dir_r is NA:
                     model_name, model_ver, mcmc_pars )
{
  if ( is.null( m ) ) { m <- model }
  if ( ! is.null( m ) ) {
    al_.file <- sprintf( "%s/al_.csv",  get_path( "model_vars_dir", m = m ) )
  } else {
    if ( is.na( model_dir_r ) ) {
      model_vars_dir <- get_path( "model_vars_dir", m = m, ps = ps, 
                                  celltype = celltype, tissue = tissue,
                                  model_name = model_name, model_ver = model_ver,
                                  mcmc_pars = mcmc_pars )
      al_.file <- sprintf( "%s/al_.csv", model_vars_dir )
      } else {
        analysis_tissue_dir <- get_path( 
          "analysis_tissue_dir", ps = ps, 
          celltype = celltype, tissue = tissue )
        al_.file <- sprintf( "%s/%s/model_vars/al_.csv", 
                             analysis_tissue_dir, model_dir_r )  
      }
    }
  
  if ( ! file.exists( al_.file ) )
    { al_ <- NULL } else { 
      al_ <- read.csv( al_.file )$x[ 1 : 6 ]
    }
  return( al_ )  
}


get_al_9 <- function( model = NULL, m = NULL,
                     ps, celltype, tissue, model_dir_r = NA,
                     # if model_dir_r is NA:
                     model_name, model_ver, mcmc_pars )
{
  if ( is.null( m ) ) { m <- model }
  if ( ! is.null( m ) ) {
    al_.file <- sprintf( "%s/al_.csv",  get_path( "model_vars_dir", m = m ) )
  } else {
    if ( is.na( model_dir_r ) ) {
      model_vars_dir <- get_path( "model_vars_dir", m = m, ps = ps, 
                                  celltype = celltype, tissue = tissue,
                                  model_name = model_name, model_ver = model_ver,
                                  mcmc_pars = mcmc_pars )
      al_.file <- sprintf( "%s/al_.csv", model_vars_dir )
    } else {
      analysis_tissue_dir <- get_path( 
        "analysis_tissue_dir", ps = ps, 
        celltype = celltype, tissue = tissue )
      al_.file <- sprintf( "%s/%s/model_vars/al_.csv", 
                           analysis_tissue_dir, model_dir_r )  
    }
  }
  
  if ( ! file.exists( al_.file ) )
  { al_9 <- NULL } else { 
    al_9 <- read.csv( al_.file )$x[ 1 : 9 ]
  }
  return( al_9 )  
}


convert_exit_flow_table_into_entry_flow_table <- function( exit.flow.table, 
                                                           al_ )
{
  if ( is.null( al_ ) | ( 
    ! ( exists( "exit.flow.table" ) && !is.null( exit.flow.table ) && 
        ncol( exit.flow.table ) == 5 ) ) ) {
    entry.flow.table <- NULL 
  } else {
    entry.flow.table <- 
      exit.flow.table %>%
      mutate( al_.entry = al_[ population.1 ] ) %>%                              
      mutate( al_.exit = al_[ population.2 ] ) %>%
      mutate( al_.entry = ifelse( al_.entry == 0, 1e-9, al_.entry ) ) %>% 
      mutate( al_.exit = ifelse( al_.exit == 0, 1e-9, al_.exit ) ) %>% 
      mutate( flow = flow * ( al_.entry / al_.exit ) ) %>%
      mutate( flow.part = "entry" ) %>% 
      dplyr::rename( population.1.old = population.1,
                     population.1 = population.2 ) %>% 
      dplyr::rename( population.2 = population.1.old ) %>% 
      dplyr::select( par, population.1, population.2, flow, flow.part ) %>% 
      arrange( population.1, population.2 )
  }
  return( entry.flow.table )
}


get_entry_flow_table <- function( ps, celltype, tissue, model_dir_r )
{
  entry.flow.table <- 
    convert_exit_flow_table_into_entry_flow_table( 
      exit.flow.table = get_exit_flow_table( 
        ps = ps, celltype = celltype, tissue = tissue, 
        model_dir_r = model_dir_r ), 
      al_ = get_al_( ps = ps, celltype = celltype, tissue = tissue,
                     model_dir_r = model_dir_r ) 
    )
  return( entry.flow.table )
}
