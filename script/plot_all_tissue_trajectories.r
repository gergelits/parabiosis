# plot_all_tissue_trajectories.r

plot_all_tissue_trajectories <- function( m_tisNA = NULL )
{
  model_id <- get_model_id( m = m_tisNA )
  fig_fname_start <- get_fig_fname_start( m = m_tisNA )
  figs_dir <- get_path( "figs_m", m = m_tisNA )
  figs_traj_dir <- get_path( "figs_m_traj", m = m_tisNA )
  
  tissues_all <- get_pmc()$tissue.all.ordered
  tissues_no_bl <- tissues_all[ -which( tissues_all == "Blood" ) ]
  tissues_bl <- c( tissues_no_bl, "Blood" )
  
  n_figs_per_A4 <- 4
  n_figs <- length( tissues_bl )
  n_A4 <- ceiling( n_figs / n_figs_per_A4 )
  ns <- list(); 
  ns$n_figs <- n_figs; ns$n_A4 <- n_A4; ns$n_figs_per_A4 <- n_figs_per_A4; 

  # modelled tissue trajectories
  l_figs_in_row_tis <- get_list_of_fig_rows( m_tisNA = m_tisNA, ns = ns, 
                                             traj_type = "tis" )
  l_trajs_multiple_pages <- plot_traj_multipage( 
    m_tisNA = m_tisNA, ns = ns, l_figs_in_row = l_figs_in_row_tis, 
    traj_type = "tis" )
  
  # other trajectories
  l_figs_in_row_pooled <- get_list_of_fig_rows( m_tisNA = m_tisNA, ns = ns, 
                                                traj_type = "pooled" )
  l_figs_in_row_blood <- get_list_of_fig_rows( m_tisNA = m_tisNA, ns = ns, 
                                               traj_type = "blood" )
  plot_traj_multipage( 
    m_tisNA = m_tisNA, ns = ns, l_figs_in_row = l_figs_in_row_pooled,
    traj_type = "pooled" )
  plot_traj_multipage(
    m_tisNA = m_tisNA, ns = ns, l_figs_in_row = l_figs_in_row_blood,
    traj_type = "blood" )
  
  if ( OUTPUT_SVG_PER_PAGES <- FALSE ) {
    cm_to_in <- 1 / 2.54
    cm_to_in2 <- cm_to_in * 2
    for ( i in 1 : n_A4 ) {
      # new version of ggpubr package needed for .svg
      # Not working on cluster?
      svg_fname <- sprintf( "%s/%s_Figure_S17_p%i.svg", 
                            figs_traj_dir, fig_fname_start, i )
      try(
        ggexport_my( 
          l_trajs_multiple_pages[[ i ]], filename = svg_fname,
          width = 21 * cm_to_in2, height = 29.7 * cm_to_in2 )
      )
    }
  }
}
