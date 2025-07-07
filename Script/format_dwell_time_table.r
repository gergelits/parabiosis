# format_dwell_time_table.r

format_dwell_time_table <- function( csv_name, ps = NA )
{
  csv_tmp <- read_load_multi( fun = "read_csv", full_path = csv_name )          
  
  d_tissues <- 
    csv_tmp %>% 
    dplyr::select( -c( dwell_time_hdi_80_lower, dwell_time_hdi_80_upper, 
                       dwell_time_hdi_50_lower, dwell_time_hdi_50_upper, 
                       par ) ) %>% 
    distinct() %>% 
    mutate( dwell_time_mode_DAYS = round( dwell_time_mode_DAYS, digits = 1 ) ) %>% 
    left_join( ., get_pmc()$dTissueAllOrderedGroup.f %>% 
                 dplyr::select( tissue.all, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) ) %>% 
    tidyr::spread( cellstate, dwell_time_mode_DAYS ) %>% 
    dplyr::rename( cd69p = `CD69+`, activ = `Antigen-experienced`, naive = Naive ) %>% 
    arrange( f.tissue.group ) %>% 
    mutate( tissue_group = f.tissue.group %>% as.character ) %>% 
    dplyr::select( celltype, tissue, tissue_group, cd69p, activ, naive, 
                   f.tissue.group )
    
    
  d_tissue_group_means <- 
    d_tissues %>% 
    group_by( f.tissue.group ) %>% 
    summarise( celltype = first( celltype ),
               tissue = "Median",
               cd69p = median( cd69p ) %>% round( digits = 1 ),
               activ = median( activ ) %>% round( digits = 1 ), 
               naive = median( naive ) %>% round( digits = 1 ) ) %>% 
    mutate( tissue_group = f.tissue.group %>% as.character ) %>% 
    dplyr::select( celltype, tissue, tissue_group, cd69p, activ, naive )
    
    
  d_tissues %>% 
    dplyr::select( celltype, tissue, tissue_group, cd69p, activ, naive ) %>% 
    rbind( ., "", "", d_tissue_group_means ) %>% 
    write_csv( file = sprintf( "%s_summary.csv", sub( '_CrI\\.csv$', '', csv_name ) ) )
}


plot_dwell_times <- function( m_tisNA = NULL, 
                              csv_name, 
                              plot = TRUE, plot_CrI = TRUE,
                              sel_models_file_name ) 
{
  celltype <- m_tisNA$celltype
  ps <- m_tisNA$ps
  mcmc_pars <- m_tisNA$mcmc_pars
  
  csv_tmp <- read_load_multi( fun = "read_csv", full_path = csv_name )
  
  dDwellTime_0 <- csv_tmp %>% 
    dplyr::select( -par ) %>% 
    distinct() %>% 
    left_join( ., get_pmc()$dTissueAllOrderedGroup.f %>% 
                 dplyr::select( tissue.all, f.tissue, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) ) %>% 
    mutate( cellstate = ifelse( 
      cellstate == "Antigen-experienced", "activated", cellstate ) ) 

  
  if ( celltype != "Treg" ) {
    dDwellTime <- dDwellTime_0 %>% 
      mutate( cellstate = ifelse( cellstate == "Naive", "naive", cellstate ) ) %>% 
      mutate( cellstate = factor( 
        cellstate, levels = c( "naive", "activated", "CD69+" ) ) )
    
  } else {
    dDwellTime <- dDwellTime_0 %>% 
      mutate( cellstate = ifelse( cellstate == "Naive", "resting", cellstate ) ) %>% 
      mutate( cellstate = factor( 
        cellstate, levels = c( "resting", "activated", "CD69+" ) ) )
  }
  
  xintercept1 <- dDwellTime$f.tissue[ 
    dDwellTime$f.tissue.group == "BoneMarrow" ] %>% unique() %>% length()
  xintercept2 <- dDwellTime$f.tissue[ 
    dDwellTime$f.tissue.group == "Lymphoid" ] %>% unique() %>% length()
  xintercept3 <- dDwellTime$f.tissue[ 
    dDwellTime$f.tissue.group == "Non-lymphoid" ] %>% unique() %>% length()
  
  cm_to_in <- 1 / 2.54   
  f <- cm_to_in * 1.7
  ylim_max <- ifelse( celltype %in% c( "Tconv" ), 270, NA )
  ylim_max <- ifelse( celltype %in% c( "B", "NK" ), 130, ylim_max )
  gg.fig.5E <- dDwellTime %>% 
    ggplot( aes( x = f.tissue, y = dwell_time_mode_DAYS, fill = cellstate ) ) +
    geom_bar( stat = "identity", position = "dodge" ) +
    geom_errorbar( aes( ymin = dwell_time_hdi_80_upper, 
                        ymax = dwell_time_hdi_80_lower ), width = 0.2,
                        position = position_dodge( 0.9 ) ) + 
    geom_vline( xintercept = c( xintercept1 + 0.5, 
                                xintercept1 + xintercept2 + 0.5,
                                xintercept1 + xintercept2 + xintercept3 + 0.5 ),
                linetype = "dashed" ) +
    labs( y = "Mean Dwell Time (days)", fill = "Cell state" ) +    
    
    ylim( 0, ylim_max ) +
  
    scale_x_discrete( labels = c( "BoneMarrow" = "Bone marrow" ) ) +
    theme_classic() +
    theme( plot.title = element_text( size = 12 * f, hjust = 0 ), 
           axis.title.x = element_blank(),
           axis.title.y = element_text( size = 12 * f, hjust = 0.5 ),  
           axis.text.x = element_text( angle = 45, hjust = 1, size = 12 * f ),
           axis.text.y = element_text( size = 12 * f ),  
           panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.text = element_text( size = 12 * f ),
           legend.title = element_text( size = 12 * f ),
           strip.text.x = element_text( size = 12 * f ) )
  
  if ( plot ) {
    fig_fname_start <- get_fig_fname_start( m = m_tisNA )

    ggexport_my( 
      gg.fig.5E, 
      filename = sprintf( "%s/%s_Figure_5E_dwell.pdf", 
                          get_path( "figs_m_dwell", m = m_tisNA ),
                          fig_fname_start ), 
      width = 21 * cm_to_in, height = 6 * cm_to_in )
    
    ggexport_my( 
      gg.fig.5E, 
      filename = sprintf( "%s/%s_Figure_5E_dwell.pdf", 
                          get_path( "figs_m", m = m_tisNA ),
                          fig_fname_start ), 
      width = 21 * cm_to_in, height = 6 * cm_to_in )
    
    # NOTE: To avoid rare error due to accessing from more threads in parallel
    try(
    add_plot_to_figs_rda( 
      m = m_tisNA,
      gg_add = gg.fig.5E, gg_add_name = "gg.fig.5E",
      gg_figs_rda_name = 
        get_figs_rda_name( fig_filename_start = fig_fname_start,
                           end = "gg_fig_5E.rda" ) )
    )
  }
}
