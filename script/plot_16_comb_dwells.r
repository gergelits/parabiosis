# plot_16_comb_dwells.r

plot_16_comb_dwells <- function( m_tisNA, n_sim = 10000, export = TRUE )
{
  l_gg_comb <- list()
  # NOTE: BoneMarrow should not be plotted here - 16 figures per A4
  tissues <- get_pmc()$dTissueAllOrderedGroup$tissue.all.ordered[ 2 : 17 ]   
  for ( tissue in tissues )
  {
    i <- which( tissue == tissues )
    l_gg_comb[[ i ]] <- ggplot() + theme_void()
    try( { 
      m <- model( m = m_tisNA, tissue = tissue )
      sim_10000 <- sim_combined_dwell( m = m, n_sim = n_sim )
      if ( ! is.null( sim_10000 ) ) {
        gg_comb <- plot_combined_dwell( m = m, sim = sim_10000, 
                                        ylim_max = 500, export = FALSE )
        l_gg_comb[[ i ]] <- gg_comb
      }
    } )
  }
  
  gg_16_comb_dwells <- ggpubr::ggarrange(
    plotlist = l_gg_comb[ 1 : 16 ], ncol = 4, nrow = 4, 
    labels = LETTERS[ 1 : length( tissues ) ] ); gg_16_comb_dwells
  
  if ( export ) {
    cm_to_in <- 1 / 2.54   
    f <- cm_to_in * 1.7
    fig_fname_start <- get_fig_fname_start( m = m_tisNA ) 
    figs_dwell_dir <- get_path( "figs_m", m = m_tisNA )
    pdf_name <- sprintf( "%s/%s_%isims_16_comb_dwell_steps_past_fut.pdf", 
                         figs_dwell_dir, fig_fname_start, n_sim )
    ggexport_my( gg_16_comb_dwells, filename = pdf_name, 
                 width = 29.7 * cm_to_in * 2, height = 21 * cm_to_in * 2 )
  }
  
  return( gg_16_comb_dwells )
}
