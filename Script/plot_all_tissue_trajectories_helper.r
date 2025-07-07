# plot_all_tissue_trajectories_helper.r

load_one_obj <- function( f )
{
  env <- new.env()
  nm <- load( f, env )[ 1 ]
  return( env[[ nm ]] )
}


get_list_of_fig_rows <- function( m_tisNA, ns, traj_type = "tis" )
{
  fig_fname_start <- get_fig_fname_start( m = m_tisNA )
  l_figs_in_row <- list()
  gg_empty_plot <- ggplot() + theme_void()
  tissues_all <- get_pmc()$tissue.all.ordered
  tissues_no_bl <- tissues_all[ -which( tissues_all == "Blood" ) ]
  tissues_bl <- c( tissues_no_bl, "Blood" )
  
  for ( i in 1 : ( ns$n_A4 * ns$n_figs_per_A4 ) ) {
    l_figs_in_row[[ i ]] <- gg_empty_plot
  }
  
  for ( i in ( 1 : ns$n_figs ) ) {
    if ( tissues_bl[ i ] != "Blood" ) {
      trajs_rda <- sprintf( "%s/%s_%s_%s_trajs.rda", get_path( 
        "model_dir", m = model( m = m_tisNA, tissue = tissues_bl[ i ] ) ), 
        fig_fname_start, tolower( tissues_bl[ i ] ), traj_type )
    } else {
      # tissues_bl[ i ] == "Blood"
      trajs_rda <- sprintf( 
        "%s/%s_%s_%s_trajs.rda", get_path( 
          "model_dir", m = model( m = m_tisNA, tissue = "Brain" ) ), 
        fig_fname_start, tolower( "Blood" ), traj_type )
    }
    
    if ( file.exists( trajs_rda ) ) {
      l_figs_in_row[[ i ]] <- read_load_multi( fun = "load_one_obj", 
                                               full_path = trajs_rda )
    } 
  }
  return( l_figs_in_row )
}


plot_traj_multipage <- function( m_tisNA, ns, l_figs_in_row, traj_type ) 
{
  fig_fname_start <- get_fig_fname_start( m = m_tisNA )
  figs_dir <- get_path( "figs_m", m = m_tisNA )
  tissues_all <- get_pmc()$tissue.all.ordered
  tissues_no_bl <- tissues_all[ -which( tissues_all == "Blood" ) ]
  tissues_bl <- c( tissues_no_bl, "Blood" )
  
  l_trajs_multiple_pages <- list()
  if ( traj_type == "tis" ) {
    labels_type <- LETTERS[ 1 : length( tissues_bl ) ]
  } else {
    labels_type <- substring( tissues_no_bl, 1, 3 )
  }
  for ( i in 1 : ns$n_A4 ) {
    js_per_i <- ( ( i - 1 ) * ns$n_figs_per_A4 ) + ( 1 : ns$n_figs_per_A4 )
    l_trajs_multiple_pages[[ i ]] <- ggpubr::ggarrange( 
      plotlist = l_figs_in_row[ js_per_i ], 
      ncol = 1, nrow = ns$n_figs_per_A4, labels = labels_type[ js_per_i ]
    )
  }
  
  cm_to_in <- 1 / 2.54
  cm_to_in2 <- cm_to_in * 2
  pdf_name <- sprintf( "%s/%s_Figure_S17_trajs_%s_multipage.pdf", 
                       figs_dir, fig_fname_start, traj_type )
  ggexport_my( l_trajs_multiple_pages, filename = pdf_name,
               width = 21 * cm_to_in2, height = 29.7 * cm_to_in2 )
  
  return( l_trajs_multiple_pages )
}
