# plot_figures_3_5.r 

plot_figures_3_5 <- function( m_tisNA = NULL )
{
  celltype <- m_tisNA$celltype; ps <- m_tisNA$ps; mcmc_pars <- m_tisNA$mcmc_pars;
  pft <- parabiosis_flow_titles( celltype = celltype )
  sel_models_file_name <- get_sel_models_file_name( m_tisNA )
  model_id <- get_model_id( sel_models_file_name )
  MIN_0 <- 0.0001
  source( sprintf( "%s/%s", ps$CODE_PATH, "plot_figures_3_5_funs.r" ), 
          local = TRUE )
  
  # Core / execution:
  dAllTissues_Models_HDI <- get_dAllTissues_Models_HDI_simple( 
    ps = ps, celltype = celltype, selected_models_file = 
      sprintf( "%s/%s", ps$ANALYSIS_PATH, sel_models_file_name ) )
  
  dAllTissues_Models_HDI_b <-
    dAllTissues_Models_HDI %>% 
    left_join( ., get_pmc()$dTissueAllOrderedGroup.f, 
               by = c( "tissue" = "tissue.all" ) ) %>%
    
    # for multi-tissue models
    mutate( par_orig = par ) %>%                                               
    mutate( par = gsub( pattern = "(q[0-9]{2,2})(_i\\[1\\])", 
                        replacement = "\\1", x = par_orig ) ) %>% 
    left_join( ., pft$dFlowsGroups, by = c( "par" = "par" ) ) %>% 
    
    filter( ! is.na( flow_between ) ) %>%  
    # do not plot (de-)differentation in blood
    filter( ! ( flow_between %in% c( "1_to_2", "2_to_3", "3_to_2" ) ) ) %>%
    filter( f.cellstate.group != "Diag" ) %>% 
    
    group_by( f.flow.subgroup, f.tissue.group.longer, f.cellstate.group ) %>%    
    mutate( mode_flow_median = median( mode ) ) %>%                              
    ungroup()
    
  
  d.pval.annot <- 
    dAllTissues_Models_HDI_b %>% 
    get.tests.table() %>% `[[`( 2 ) %>%
    get.pval.annot.table()
  
  cm_to_in <- 1 / 2.54   
  f <- cm_to_in * 1.7
  
  ggfig_3_5 <- 
    dAllTissues_Models_HDI_b %>% 
    mutate( mode = ifelse( mode == 0, MIN_0, mode ) ) %>%
    mutate( mode_flow_median = 
              ifelse( mode_flow_median == 0, MIN_0, mode_flow_median ) ) %>% 
    
    ggplot( aes( x = f.tissue.group.longer, y = mode * 1000 ) ) +
    geom_jitter( aes( shape = f.cellstate.group, colour = f.cellstate.group ),
                 position = position_dodge( 0.3 ), size = 3 * f ) +
    geom_errorbar( aes( ymax = mode_flow_median * 1000, 
                        ymin = mode_flow_median * 1000,
                        colour = f.cellstate.group, 
                        width = 0.8 ) ) +
    
    scale_y_log10( n.breaks = 8,                                                 
                   labels = scales::trans_format( 
                     "log10", scales::math_format( 10^.x ) ) ) +
    annotation_logticks( sides = "l" ) +
    ylab( "Rate (events/1000 cells/day)" ) +     
  
    facet_grid( . ~ f.flow.subgroup ) +
    geom_text( data = d.pval.annot, 
               aes( label = w.p.value.text ), size = 4 * f ) +
    
    theme_classic() +
    theme( plot.title = element_text( hjust = 0 ),
           strip.text.x = element_text( size = 12 * f ),
           axis.title.x = element_blank(),
           axis.title.y = element_text( size = 12 * f ), 
           axis.text.x = element_text( angle = 45, hjust = 1, size = 12 * f ),     
           axis.text.y = element_text( size = 12 * f ),  
           legend.position = c( 0.90, 0.85 ),
           legend.title = element_blank(),
           legend.text = element_text( size = 12 * f ),
           panel.border = element_blank(),                
           panel.grid.major = element_blank(),            
           panel.grid.minor = element_blank() ); ggfig_3_5
  
  figs_dir <- get_path( "figs_m", m = m_tisNA )
  if ( !dir.exists( figs_dir ) ) { dir.create( figs_dir, recursive = TRUE ) }
  fig_filename_start <- get_fig_fname_start( m = m_tisNA )
  pdf_name <- sprintf( "%s/%s_Figure_3or5_SideBySide.pdf",
                       figs_dir, fig_filename_start )
  ggexport_my( ggfig_3_5, filename = pdf_name,
               width = 21 * cm_to_in, height = 12 * cm_to_in )
}
