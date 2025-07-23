# plot_figures_2_4.r

get_i_flow_ggfig <- function( dTMP.2, i_flow, MIN_0 = 1e-6 )
{
  # vertical line separaterors
  xintercept1 <- dTMP.2$f.tissue[ dTMP.2$f.tissue.group == "BoneMarrow" ] %>% 
    unique() %>% length()
  xintercept2 <- dTMP.2$f.tissue[ dTMP.2$f.tissue.group == "Lymphoid" ] %>% 
    unique() %>% length()
  xintercept3 <- dTMP.2$f.tissue[ dTMP.2$f.tissue.group == "Non-lymphoid" ] %>% 
    unique() %>% length()
  celltype <- dTMP.2$celltype %>% unique()
  
  hdi.80.lower.flow.subgroup.min.this <- 
    dTMP.2 %>% filter( flow_between == i_flow ) %>% 
    dplyr::select( hdi.80.lower.flow.subgroup.min ) %>% min()
  if ( hdi.80.lower.flow.subgroup.min.this == 0 ) 
    hdi.80.lower.flow.subgroup.min.this <- MIN_0
  
  hdi.80.upper.flow.subgroup.max.this <- 
    dTMP.2 %>% filter( flow_between == i_flow ) %>% 
    dplyr::select( hdi.80.upper.flow.subgroup.max ) %>% max()
  if ( hdi.80.upper.flow.subgroup.max.this == 0 ) 
    hdi.80.upper.flow.subgroup.max.this <- MIN_0
  
  cm_to_in <- 1 / 2.54   
  f <- cm_to_in * 1.7
  
  ggfig <- dTMP.2 %>% 
    filter( !is.na( f.tissue ) ) %>% 
    filter( flow_between == i_flow ) %>%
    mutate(
      hdi.80.upper = ifelse( hdi.80.upper < MIN_0, MIN_0, hdi.80.upper ),
      hdi.50.upper = ifelse( hdi.50.upper < MIN_0, MIN_0, hdi.50.upper ),
      mode = ifelse( mode < MIN_0, MIN_0, mode ),
      hdi.80.lower = ifelse( hdi.80.lower < MIN_0, MIN_0, hdi.80.lower ),
      hdi.50.lower = ifelse( hdi.50.lower < MIN_0, MIN_0, hdi.50.lower ) ) %>% 
    
    ggplot( aes( x = f.tissue, y = mode * 1e3 ) ) +                              
    geom_point( size = 3, alpha = 0.5 ) +
    scale_y_log10( n.breaks = 6,
                   limits = c( hdi.80.lower.flow.subgroup.min.this * 1e3 ,  
                               hdi.80.upper.flow.subgroup.max.this * 1e3 ),    
                   labels = scales::trans_format( 
                     "log10", scales::math_format( 10^.x ) ) ) +               
    annotation_logticks( sides = "l" ) +
    geom_crossbar( aes( ymin = hdi.80.lower * 1e3, ymax = hdi.80.upper * 1e3 ), 
                   colour = "grey", width = 0.6 ) +
    geom_crossbar( aes( ymin = hdi.50.lower * 1e3, ymax = hdi.50.upper * 1e3 ), 
                   colour = "blue" ) +
    
    # graphics
    geom_vline( xintercept = 
                  c( xintercept1 + 0.5, 
                     xintercept1 + xintercept2 + 0.5,
                     xintercept1 + xintercept2 + xintercept3 + 0.5 ) ) +
    geom_hline( yintercept = 1, colour = "darkred", 
                linetype = 2, linewidth = 1 )  +
    labs( title = dTMP.2$flow.title[ dTMP.2$flow_between == i_flow ] %>% 
            na.omit %>% unique(),
          y = "Rate (events/1000 cells/day)" ) +     
    scale_x_discrete( labels = c( "BoneMarrow" = "Bone marrow" ) ) +
    theme_classic() +
    theme( plot.title = element_text( size = 12 * f, hjust = 0 ), 
           axis.title.x = element_blank(),
           axis.title.y = element_text( size = 12 * f, hjust = 1 ),  
           
           axis.text.x = element_text( angle = 45, hjust = 1, size = 12 * f ),
           axis.text.y = element_text( size = 12 * f ),  
           
           legend.text = element_text( size = 12 * f ),
           legend.title = element_text( size = 12 * f ),
           
           panel.border = element_blank(),             
           panel.grid.major = element_blank(),         
           panel.grid.minor = element_blank() ); ggfig
  
  return( ggfig )
}


plot_figures_2_4 <- function( m_tisNA = NULL, qij_I = 1L )
{
  celltype <- m_tisNA$celltype; ps <- m_tisNA$ps; mcmc_pars <- m_tisNA$mcmc_pars;
  pft <- parabiosis_flow_titles( celltype = celltype )
  sel_models_file_name <- get_sel_models_file_name( m_tisNA )
  dAllTissues_Models_HDI <- get_dAllTissues_Models_HDI_simple( 
    ps = ps, celltype = celltype, selected_models_file = 
      sprintf( "%s/%s", ps$ANALYSIS_PATH, sel_models_file_name ) )               
  
  model_id <- get_model_id( sel_models_file_name )
  
  # which part of model: 1...tissue, 2...pooled others
  gsub_pattern <- sprintf( "(q[0-9]{2,2})(_i\\[%i\\])", qij_I )
  
  dTMP.2 <- 
    dAllTissues_Models_HDI %>% 
    # For the vertical lines separating tissue groups:
    left_join( get_pmc()$dTissueAllOrderedGroup.f %>%                                   
                 dplyr::select( tissue.all, f.tissue, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) ) %>% 
    # multi-tissue models
    mutate( par_orig = par ) %>%                                                 
    mutate( par = gsub( pattern = gsub_pattern, 
                        replacement = "\\1", x = par_orig ) ) %>% 
    
    left_join( ., pft$dFlowsGroups %>%
                 dplyr::select( par, f.flow.subgroup, flow_between ),
               by = c( "par" = "par" ) ) %>%
    # get titles   
    left_join( ., pft$dFlowsTitles %>% dplyr::select( flow_between, flow.title ),
               by = c( "flow_between" = "flow_between" ) ) %>%
    
    # f.flow.subgroup - differentiation, entry, dediff, 
    mutate( f.flow.subgroup.to.plot = ifelse( 
      flow_between %in% c( "1_to_2", "4_to_1", "5_to_2", "6_to_5" ), 
      "not.to.plot", f.flow.subgroup ) ) %>% 
    group_by( f.flow.subgroup.to.plot ) %>%                                      
    mutate( 
      hdi.80.lower.flow.subgroup.min = min( hdi.80.lower, na.rm = TRUE ),
      hdi.80.upper.flow.subgroup.max = max( hdi.80.upper, na.rm = TRUE ) ) %>%
    ungroup()
    
  l_ggfigs <- list()
  FLOWS <- dTMP.2$flow_between %>% na.omit %>% unique %>% sort
  for ( i_flow in FLOWS ) { 
    ggfig <- get_i_flow_ggfig( dTMP.2 = dTMP.2, i_flow = i_flow )
    l_ggfigs[[ which( i_flow == FLOWS ) ]] <- ggfig
  }
  
  fig_filename_start <- get_fig_fname_start( m = m_tisNA )
  figs_dir <- get_path( "figs_m", m = m_tisNA, ps = ps, celltype = celltype )
  if ( ! dir.exists( figs_dir ) ) { dir.create( figs_dir, recursive = TRUE ) }
  
  cm_to_in <- 1 / 2.54   
  f <- cm_to_in * 1.7
  if ( ( FIG_5BCD <- TRUE ) & ( qij_I == 1 ) ) {
    gg.fig.5BCD <- ggpubr::ggarrange( 
      l_ggfigs[[ which( "1_to_4" == FLOWS ) ]], 
      l_ggfigs[[ which( "2_to_5" == FLOWS ) ]], 
      l_ggfigs[[ which( "3_to_6" == FLOWS ) ]],
      labels = c( "B", "C", "D" ),
      ncol = 3, nrow = 1 )
    
    pdf_name <- sprintf( "%s/%s_Figure_5BCD.pdf", figs_dir, fig_filename_start )
    ggexport_my( gg.fig.5BCD, filename = pdf_name,
                 width = 21 * cm_to_in, height = 6 * cm_to_in )
    
    # SOMETIMES due to accessing from multiple threads in parallel?
    # Error in load(full_gg_figs_rda_name) : error reading from connection
    try( 
      add_plot_to_figs_rda( 
        m = m_tisNA,
        gg_add = gg.fig.5BCD, gg_add_name = "gg.fig.5BCD",
        gg_figs_rda_name = 
          get_figs_rda_name( fig_filename_start = fig_filename_start,
                             end = "gg_fig_5BCD.rda" ) )
    )
  }
  
  if ( FIG_DIAG_18_PARAMS <- TRUE ) {
    gg.fig.one.A4.debug4 <- ggpubr::ggarrange( 
      l_ggfigs[[ which( "1_to_4" == FLOWS ) ]], 
      l_ggfigs[[ which( "2_to_5" == FLOWS ) ]], 
      l_ggfigs[[ which( "1_to_2" == FLOWS ) ]], 
      l_ggfigs[[ which( "2_to_3" == FLOWS ) ]],
      
      l_ggfigs[[ which( "4_to_5" == FLOWS ) ]], 
      l_ggfigs[[ which( "5_to_6" == FLOWS ) ]], 
      l_ggfigs[[ which( "4_to_1" == FLOWS ) ]], 
      l_ggfigs[[ which( "3_to_2" == FLOWS ) ]],
      
      l_ggfigs[[ which( "5_to_1" == FLOWS ) ]], 
      l_ggfigs[[ which( "6_to_1" == FLOWS ) ]], 
      l_ggfigs[[ which( "5_to_2" == FLOWS ) ]], 
      l_ggfigs[[ which( "3_to_6" == FLOWS ) ]],
      
      l_ggfigs[[ which( "q17" == FLOWS )    ]], 
      l_ggfigs[[ which( "q28" == FLOWS )    ]], 
      l_ggfigs[[ which( "q39" == FLOWS )    ]], 
      l_ggfigs[[ which( "brf_free_param" == FLOWS ) ]], 
      
      l_ggfigs[[ which( "6_to_3" == FLOWS ) ]], 
      l_ggfigs[[ which( "6_to_5" == FLOWS ) ]], 
      NULL, 
      NULL,
      labels = c( "A", "B", "q12", "q23", 
                  "C", "D", "q41", "q32",
                  "E", "F", "q52", "q36",
                  "hdr1", "hdr2", "hdr3", "brf_free_param",
                  "q63", "q65" ),
      ncol = 4, nrow = 5 )
    
    pdf_name <- sprintf( "%s/%s_diag_18_params_part%i.pdf", 
                         figs_dir, fig_filename_start, qij_I ) 
    
    ggexport_my( gg.fig.one.A4.debug4, filename = pdf_name, 
                 width = 16, height = 15 )
  }
}  



