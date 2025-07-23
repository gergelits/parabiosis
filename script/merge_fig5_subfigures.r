# merge_fig5_subfigures.r 

merge_fig5_subfigures <- function( m_tisNA = NULL )
{
  celltype <- m_tisNA$celltype; ps <- m_tisNA$ps; mcmc_pars <- m_tisNA$mcmc_pars;
  figs_dir <- get_path( "figs_m", m = m_tisNA )
  figs_rda_dir <- get_path( "figs_m_rda", m = m_tisNA )
  fig_filename_start <- get_fig_fname_start( m = m_tisNA )           
  
  # SOMETIMES due to accessing from multiple threads in parallel?
  # Error in load(full_gg_figs_rda_name) : error reading from connection
  # NOTE: sometimes .rds / .rda was corrupted (due to cloud storing?)
  cm_to_in <- 1 / 2.54
  try({
    gg_5A_rda <- sprintf( "%s/%s", figs_rda_dir, get_figs_rda_name( 
      fig_filename_start = fig_filename_start, end = "gg_fig_5A.rda" ) )
    load( gg_5A_rda )  
    gg.fig.5A <- gg_figs_all$gg.fig.5A
  
    gg_5BCD_rda <- sprintf( "%s/%s", figs_rda_dir, get_figs_rda_name( 
      fig_filename_start = fig_filename_start, end = "gg_fig_5BCD.rda" ) )
    load( gg_5BCD_rda )
    gg.fig.5BCD <- gg_figs_all$gg.fig.5BCD
  
    gg_5E_rda <- sprintf( "%s/%s", figs_rda_dir, get_figs_rda_name( 
      fig_filename_start = fig_filename_start, end = "gg_fig_5E.rda" ) )
    load( gg_5E_rda )
    gg.fig.5E <- gg_figs_all$gg.fig.5E
    
    gg_figs.ABCDE <- ggpubr::ggarrange(
      gg.fig.5A, gg.fig.5BCD, gg.fig.5E, 
      labels = c( "A", "", "E" ), ncol = 1, nrow = 3 )
    
    pdf_name <- sprintf( "%s/%s_Figure_5ABCDE.pdf", figs_dir, 
                         fig_filename_start )
    ggexport_my( gg_figs.ABCDE, filename = pdf_name,
                 width = 21 * cm_to_in, height = 18 * cm_to_in )
  })
  
  # new ggpubr package needed for .svg (on cluster)
  try({
    svg_name <- sprintf( "%s/%s_Figure_5ABCDE.svg", 
                         figs_dir, fig_filename_start )
    ggexport_my( gg_figs.ABCDE, filename = svg_name,
                 width = 21 * cm_to_in, height = 18 * cm_to_in )
    }) 
}
