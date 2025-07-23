# prep_data_figures_6.r                               

write_aggregated_group_flows <- function( m_tisNA = NULL, write_csv = FALSE )
{
  celltype <- m_tisNA$celltype; ps <- m_tisNA$ps; mcmc_pars <- m_tisNA$mcmc_pars;
  model_id <- get_model_id( m = m_tisNA )
  model_name <- m_tisNA$model_name; model_ver <- m_tisNA$model_ver
  
  avg_flows_exit <- get_tissue_average_flows( 
    m_tisNA = m_tisNA, ps = ps, celltype = celltype, flow.part = "exit",
    model_name = model_name, model_ver = model_ver, 
    mcmc_pars = mcmc_pars ) 
  avg_flows_entry <- get_tissue_average_flows( 
    m_tisNA = m_tisNA, ps = ps, celltype = celltype, flow.part = "entry",
    model_name = model_name, model_ver = model_ver, 
    mcmc_pars = mcmc_pars )
  average.flows.celltype <- bind_rows( avg_flows_exit, avg_flows_entry )
    
  
  if ( all( dim( average.flows.celltype ) == c( 0, 0 ) ) ) { 
    return( NULL )
  }
  
  aggregated.group.flows <- 
    average.flows.celltype %>%
    mutate( arrow.group = sprintf( 
      "g_%s_%s", 
      ifelse( population.1 < population.2, population.1, population.2 ),
      ifelse( population.1 > population.2, population.1, population.2 ) ) ) %>% 
    
    mutate( median.flow.1000.round = round( 1000 * median.flow, 0 ) ) %>% 
    dplyr::select( f.tissue.group, arrow.group, population.1, population.2, par, 
                   mean.flow, flow.mean.w.by.tissue.in.tissue.group, flow.part, 
                   median.flow, median.flow.1000.round, 
                   log10.median.flow, log10.median.flow.1000 ) %>% 
    
    mutate( log10.median.flow.1000.tmp.3 = 
              ifelse( log10.median.flow.1000 < 1, 1, 
                      ifelse( log10.median.flow.1000 > 4, 
                              4, log10.median.flow.1000 ) ) - 1 ) %>% 
    mutate( log10.median.flow.1000.deg.3 = 
              log10.median.flow.1000.tmp.3 / 3 * 90 ) %>% 
    
    mutate( log10.median.flow.1000.tmp.4 = 
              ifelse( log10.median.flow.1000 < 0, 0, 
                      ifelse( log10.median.flow.1000 > 4, 
                              4, log10.median.flow.1000 ) ) ) %>% 
    mutate( log10.median.flow.1000.deg.4 = 
              log10.median.flow.1000.tmp.4 / 4 * 90 ) %>% 
    
    dplyr::select( -c( 
      log10.median.flow.1000.tmp.3, log10.median.flow.1000.tmp.4,
      log10.median.flow.1000.deg.4 ) ) %>% 
    
    arrange( arrow.group, f.tissue.group, par )
    
  
  if ( write_csv ) {
    fig_filename_start <- get_fig_fname_start( m = m_tisNA )
    flow_fig_dir <- get_path( "flow_fig", ps = ps, celltype = celltype )
    if( !dir.exists( flow_fig_dir ) ) { dir.create( flow_fig_dir ) }
    
    write_csv( 
      aggregated.group.flows, 
      sprintf( "%s/%s_flow_diagram_medians_means_degs.csv", 
               flow_fig_dir, fig_filename_start ) )
    
    write_csv(
      aggregated.group.flows %>% filter( f.tissue.group != "BoneMarrow" ), 
      sprintf( "%s/%s_flow_diagram_medians_means_degs_noBM.csv", 
               flow_fig_dir, fig_filename_start ) )
    }
  return( aggregated.group.flows )
}


get_cellstate_areas <- function( 
    ps = PS, f_tissue = "Brain", celltype = "Treg", 
    f_week = 0, hd  = "host" ) 
{
  stopifnot( hd %in% c( "host", "donor" ) )
  
  parabio.file.name <- 
    sprintf( "%s/%s/%s_parabiosis_data_%s.csv", 
    ps$PROCESSED_PATH, celltype, celltype, 
    tolower( ifelse( f_tissue == "Blood", "Brain", f_tissue ) ) )
  parabio.data <- read.csv( parabio.file.name )
  colnames( parabio.data ) <- 
    gsub( pattern = "\\.(.*)\\.", 
          replacement = ".celltype.", colnames( parabio.data ) )
  if ( hd == "host" ) {
    parabio.data$sel.hd.celltype.naive <- parabio.data$host.celltype.naive
    parabio.data$sel.hd.celltype.activ <- parabio.data$host.celltype.activ
    parabio.data$sel.hd.celltype.cd69p <- parabio.data$host.celltype.cd69p
  } else {
    parabio.data$sel.hd.celltype.naive <- parabio.data$donor.celltype.naive
    parabio.data$sel.hd.celltype.activ <- parabio.data$donor.celltype.activ
    parabio.data$sel.hd.celltype.cd69p <- parabio.data$donor.celltype.cd69p
  }
    
  areas_cellstates <-
    parabio.data %>%
    filter( tissue == f_tissue, week == f_week ) %>% 
    summarise( area_naive = mean( sel.hd.celltype.naive ), 
               area_activ = mean( sel.hd.celltype.activ ), 
               area_cd69p = mean( sel.hd.celltype.cd69p ) ) %>%  
    sweep( ., 1, rowSums( . ), "/" )
  
  dCellstateAreas <- tibble( tissue = f_tissue, week = f_week, hd = hd,
          cellstate = sub( "area_", "", names( areas_cellstates ) ),
          area = as.numeric( areas_cellstates ),
          radius = sqrt( area ) )

  return( dCellstateAreas )
}


get_aggr_tissuegroup_cellstate_areas <- function( 
    m_tisNA = NULL, write_csv = FALSE )
{
  ps <- m_tisNA$ps; celltype <- m_tisNA$celltype; 
  dTissueAllOrderedGroup.f <- get_pmc()$dTissueAllOrderedGroup.f
  
  d.tissue.all <- NULL
  for( tissue in dTissueAllOrderedGroup.f$tissue.all ) {
    d.tissue.all <- d.tissue.all %>% 
      bind_rows( ., get_cellstate_areas( ps, tissue, celltype ) )
  }
  
  median_tissuegroup_cellstate_areas <- 
    d.tissue.all %>% 
    left_join( ., dTissueAllOrderedGroup.f %>% 
                 dplyr::select( tissue.all, tissue.group ), 
               by = c( "tissue" = "tissue.all") ) %>% 
    group_by( tissue.group, cellstate ) %>% 
    summarise( median_tissue_group_area = median( area ), .groups = "drop" ) %>% 
    mutate( radius = sqrt( median_tissue_group_area ) ) %>% 
    # add Blood areas
    bind_rows( ., get_cellstate_areas( ps, "Blood", celltype ) %>%
                 dplyr::rename( median_tissue_group_area = area,
                                tissue.group = tissue ) ) %>% 
    dplyr::select( tissue.group, cellstate, 
                   median_tissue_group_area, radius ) %>% 
    arrange( tissue.group, cellstate )
  
  if ( write_csv ) { 
    flow_fig_dir <- get_path( "flow_fig", ps = ps, celltype = celltype )
    if ( ! dir.exists( flow_fig_dir ) ) { dir.create( flow_fig_dir ) }
    write_csv( median_tissuegroup_cellstate_areas, sprintf( 
      "%s/%s_tissuegroup_cellstates_areas.csv", flow_fig_dir, celltype ) )
  }
  
  return( median_tissuegroup_cellstate_areas )
}


get_aggr_total_counts_areas <- function( m_tisNA = NULL, 
                                         write_csv = FALSE, 
                                         count_type = "Median" )
{
  ps <- m_tisNA$ps; celltype <- m_tisNA$celltype; 
  total_counts <- read_csv( 
    sprintf( "%s/Total_counts/parabiosis_model_input_%s_counts.csv", 
             ps$PROCESSED_PATH, celltype ) )
  
  dTotalCountsAreas <-
    total_counts %>% 
    filter( Tissue %in% get_pmc()$tissue.all.ordered ) %>% 
    left_join( ., get_pmc()$dTissueAllOrderedGroupBlood.f %>% 
                 dplyr::select( tissue.all, tissue.group ), 
               by = c( "Tissue" = "tissue.all" ) ) %>% 
    # the following grouping needed for correct assignment to Mean_or_Median
    group_by( Tissue ) %>% 
    mutate( Mean_or_Median = 
              ifelse( count_type == "Median", Median, Mean ) ) %>% 
    group_by( tissue.group ) %>% 
    summarise( median_total_count = median( Mean_or_Median ),
               mean_total_count = mean( Mean_or_Median ) ) %>% 
    mutate( max_median_total_count = max( median_total_count ),
            rel_to_max_total_count = 
              median_total_count / max_median_total_count,
            radius_rel_total_count = sqrt( rel_to_max_total_count ) )
  
  if ( write_csv ) {   
    flow_fig_dir <- get_path( "flow_fig", m = m_tisNA )
    if ( !dir.exists( flow_fig_dir ) ) { dir.create( flow_fig_dir ) }
    write_csv( dTotalCountsAreas, sprintf( 
      "%s/%s_total_count_areas.csv", flow_fig_dir, celltype ) )
  }

  return( dTotalCountsAreas )
}


plot_aggr_tissuegroup_cellstate_barplot <- function( 
    m_tisNA = NULL, plot = TRUE )
{
  celltype <- m_tisNA$celltype; ps <- m_tisNA$ps
  dCellstatesProps <- NULL
  for ( tissue in get_pmc()$tissue.all.ordered )
    for ( hd in c( "host", "donor" ) )
      for ( week in c( 0, 1, 2, 4, 8, 12 ) ) {
        dCellstatesProps <- 
          dCellstatesProps %>% 
          bind_rows( ., get_cellstate_areas( 
            ps = PS, f_tissue = tissue, celltype = celltype, 
            hd = hd, f_week = week ) ) 
      }
  
  dCellstatesPropsGroups_0 <-
    dCellstatesProps %>% 
    filter( ! ( week == 0 & hd == "donor" ) ) %>%
    filter( tissue != "BoneMarrow" ) %>% 
    left_join( ., get_pmc()$dTissueAllOrderedGroupBlood.f,
               by = c( "tissue" = "tissue.all" ) ) %>% 
    group_by( tissue.group, hd, week, cellstate ) %>% 
    summarise( mean_prop = mean( area ), .groups = "drop_last" ) %>%
    ungroup() %>% 
    # needs to keep celltype
    bind_rows( ., tibble( 
      mean_prop = 1, hd = "donor", week = 0, 
      tissue.group = c( "Blood", "Non-lymphoid", "Lymphoid", "GALT" ) ) ) %>% 
    mutate( hd_wk = sprintf( "%s_wk%s", hd, str_pad( week, 2, pad = "0" ) ) ) %>%
    mutate( hd_ordered = ifelse( hd == "host", "01_host", "02_donor" ) ) %>%
    arrange( hd_ordered, week ) %>% 
    
    mutate( hd_hack = ifelse( hd == "host", "", " " ) ) %>% 
    mutate( hd_wk_ordered = factor( sprintf(
      "%s_wk%s", hd_ordered, str_pad( week, 2, pad = "0" ) ),
      labels = unique( sprintf( "%s%s%s", hd_hack, week, hd_hack ) ) ) ) %>%
  
    mutate( cellstate = ifelse( 
      cellstate == "naive" & celltype == "Treg", "resting", cellstate ) ) %>% 
    mutate( cellstate = ifelse( 
      cellstate == "activ", "activated", cellstate ) ) %>% 
    mutate( cellstate = ifelse( 
      cellstate == "cd69p", "CD69+", cellstate ) ) %>% 
    mutate( tissue.group = factor( 
      tissue.group, levels = c( "Blood", "Lymphoid", "Non-lymphoid", "GALT" ) ) )
  
  
  if ( celltype != "Treg" ) {
    dCellstatesPropsGroups <- dCellstatesPropsGroups_0 %>% 
      mutate( cellstate = factor( 
        cellstate, levels = c( "naive", "activated", "CD69+" ) ) ) 
  } else {
    dCellstatesPropsGroups <- dCellstatesPropsGroups_0 %>% 
      mutate( cellstate = factor( 
        cellstate, levels = c( "resting", "activated", "CD69+" ) ) ) 
  }
  
  cm_to_in <- 1 / 2.54   
  f <- cm_to_in * 1.7
  
  ggfig <- 
    dCellstatesPropsGroups %>% 
    ggplot( ., aes( x = hd_wk_ordered, y = mean_prop, fill = cellstate ) ) + 
    geom_bar( position = "fill", stat = "identity" ) +
    facet_grid( cols = vars( tissue.group ) ) +
    theme_classic() +
    theme( plot.title = element_text( hjust = 0 ), 
           axis.title.x = element_blank(),
           axis.title.y = element_text( size = 12 * f, hjust = 0.5 ),  
           
           axis.text.x = element_text( size = 12 * f ),
           axis.text.y = element_text( size = 12 * f ),  
           
           panel.border = element_blank(),             
           panel.grid.major = element_blank(),         
           panel.grid.minor = element_blank(),
      
           legend.text = element_text( size = 12 * f ),
           legend.title = element_text( size = 12 * f ),
           strip.text.x = element_text( size = 12 * f, vjust = 0 ),
           strip.background = element_blank(),
           plot.tag = element_text( size = 12 * f, hjust = 0 ),
           plot.tag.position = c( 0.00, 0.037 ),
           plot.margin = margin( 5.5, 5.5, 11.5, 5.5, unit = "pt" ) ) +   
      
    scale_fill_discrete( na.value = "grey" ) + 
    
    annotate( "text", x = 3.5, y = -0.20, label = "host", size = 4 * f ) +
    annotate( "text", x = 9.5, y = -0.20, label = "donor", size = 4 * f ) +
    coord_cartesian( ylim = c( 0, NA ), clip = "off" ) +
    labs( y = "Cell state proportion", fill = "Cell state", tag = "week" )

  
  if ( plot ) {
    fig_filename_start <- get_fig_fname_start( m = m_tisNA )
    pdf_name <- sprintf( 
      "%s/%s_data_Figure_5A_cellstates_barchart.pdf", 
      get_path( "figs", ps = ps, celltype = celltype ), celltype )
    ggexport_my( ggfig, filename = pdf_name, 
                 width = 21 * cm_to_in, height = 6 * cm_to_in )
    
    ## SOMETIMES due to accessing from multiple threads in parallel?
    # Error in load(full_gg_figs_rda_name) : error reading from connection
    try(
      add_plot_to_figs_rda( 
        m = m_tisNA,
        gg_add = ggfig, gg_add_name = "gg.fig.5A",
        gg_figs_rda_name = get_figs_rda_name( 
          fig_filename_start = fig_filename_start, end = "gg_fig_5A.rda" ) )
    )
  }
  return( ggfig )
}
