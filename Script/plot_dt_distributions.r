# plot_dt_distributions.r

get_dwell_times_table <- function( 
    m_tisNA = NULL, celltype, ps, model_id )
{ 
  csv_full_path <- sprintf( 
    "%s/%s_dwell_summary.csv",
    get_path( "figs_m", m = m_tisNA ),
    get_fig_fname_start( m = m_tisNA ) )

  if( !file.exists( csv_full_path ) ) { return( NULL ) }
  
  d_dwell_times_tmp <- read_load_multi( fun = "read.csv", 
                                        full_path = csv_full_path )
  d_dwell_times <- 
    d_dwell_times_tmp %>% 
    filter( tissue %in% get_pmc()$dTissueAllOrderedGroup$tissue.all.ordered )
    
  return( d_dwell_times )
}


plot_dwell_time_density <- function( 
    mean_dt, probs = c( 0.01, 0.25, 0.5, 0.75, 0.90 ), 
    tissue, all_tissues_max = 0 )
{
  lambda <- 1 / mean_dt 
  qs <- qexp( p = probs, rate = lambda )
  q975 <- qexp( p = 0.975, rate = lambda )
  all_tissues_max <- max( all_tissues_max, q975 * 1.1 )
  
  cm_to_in <- 1 / 2.54   
  f <- cm_to_in * 1.7
  p <- ggplot(
    data.frame( 
      x = c( qs, all_tissues_max ),
      q_label = c( sprintf( "q%s=%s", probs * 100, round( qs ) ), "" ) ), 
    aes( x ) ) + 
    ylim( 0, 0.1 ) + 
    xlim( 0, 150 ) +
    stat_function( fun = function( x ) dexp( x, rate = lambda ) ) +
    labs( title = tissue, 
          x = "Dwell time (days)", y = "Probability density" ) + 
    geom_text( aes( y = order( -x ) / 200, label = q_label ), 
               colour = "red", size = 4 * f ) +
    geom_vline( xintercept = qs ) +
    geom_text( y = 0, x = mean_dt, label = sprintf( "mean=%s", mean_dt ), 
               colour = "darkred", size = 4 * f ) +
    geom_vline( xintercept = mean_dt, linetype = "dashed" ) +
    theme( plot.title = element_text( size = 12 * f ),
           strip.text.x = element_text( size = 12 * f ),
           axis.title.x = element_text( size = 12 * f ),
           axis.title.y = element_text( size = 12 * f ), 
           axis.text.x = element_text( size = 12 * f ),     
           axis.text.y = element_text( size = 12 * f ) )
  return( p )
}


plot_dt_distributions <- function( m_tisNA = NULL, celltype, cellstate, ps,
                                   model_id, mcmc_pars ) 
{
  d_dwell_times <- get_dwell_times_table( 
    m_tisNA = m_tisNA,
    celltype = m_tisNA$celltype, ps = m_tisNA$ps, 
    model_id = m_tisNA$model_id )
  
  l_ggfigs <- list()
  n_tissues <- nrow( d_dwell_times )
  for ( i in 1 : n_tissues ) {
    l_ggfigs[[ i ]] <- ( plot_dwell_time_density( 
      mean = d_dwell_times[[ cellstate ]][ i ], 
      tissue = d_dwell_times$tissue[ i ],
      all_tissues_max = 150 ) )
  }
  gg_empty_plot <- ggplot() + theme_void()
  for ( i in ( n_tissues + 1 ) : 18 ) {
    l_ggfigs[[ i ]] <- gg_empty_plot
  }
  
  n_fig_per_A4 <- 3 * 3
  gg.fig.one.A4.1 <- ggpubr::ggarrange(
    plotlist = l_ggfigs[ 1 : n_fig_per_A4 ],
    labels = c( LETTERS[ 1 : n_fig_per_A4 ] ),
    ncol = 3, nrow = 3 )
  
  gg.fig.one.A4.2 <- ggpubr::ggarrange(
    plotlist = l_ggfigs[ ( n_fig_per_A4 + 1 ) : ( 2 * n_fig_per_A4 ) ],
    labels = c( LETTERS[ ( n_fig_per_A4 + 1 ) : n_tissues ] ),
    ncol = 3, nrow = 3 )
  
  l_dt_per_state <- list()
  l_dt_per_state[[ 1 ]] <- gg.fig.one.A4.1
  l_dt_per_state[[ 2 ]] <- gg.fig.one.A4.2
  
  figs_dir <- get_path( "figs_m_dwell", m = m_tisNA )
  fig_filename_start <- get_fig_fname_start( m = m_tisNA  )
  
  cm_to_in <- 1 / 2.54   
  f <- cm_to_in * 1.7
  pdf_name <-  sprintf( "%s/%s", figs_dir, sprintf( 
    "%s_Figure_S18_dwell_%s.pdf", fig_filename_start, cellstate ) )
  ggexport_my( l_dt_per_state, filename = pdf_name,
               width = 21 * cm_to_in * 1.5, height = 18 * cm_to_in * 1.5 )
}
