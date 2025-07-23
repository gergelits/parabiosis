# plot_figure_combined_4.r 

get_fig_3_5_rates_all <- function( m_tisNA = NULL )
{
  source( sprintf( "%s/%s", m_tisNA$ps$CODE_PATH, "plot_figures_3_5_funs.r" ), 
          local = TRUE )
  
  celltypes <- c( "Tconv", "CD8", "Treg" )
  dAllTissues_Models_HDI_b_all <- NULL
  for ( celltype in celltypes ) {
    m_tisNA$celltype <- celltype
    ps <- m_tisNA$ps; mcmc_pars <- m_tisNA$mcmc_pars;
    pft <- parabiosis_flow_titles( celltype = celltype )
    sel_models_file_name <- get_sel_models_file_name( m_tisNA )
    model_id <- get_model_id( sel_models_file_name )
    selected_models_file <- sprintf( 
      "%s/%s", ps$ANALYSIS_PATH, sel_models_file_name )
    
    if ( file.exists( selected_models_file ) ) {
      dAllTissues_Models_HDI <- get_dAllTissues_Models_HDI_simple( 
        ps = ps, celltype = celltype, 
        selected_models_file = selected_models_file )
      
      dAllTissues_Models_HDI_b <- dAllTissues_Models_HDI %>% 
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
        filter( f.flow.subgroup %in% c( 
          "Apoptosis", "Tissue entry\n(rate rel. to blood counts)" ) ) %>% 
        group_by( f.flow.subgroup, f.tissue.group.longer, f.cellstate.group ) %>%    
        mutate( mode_flow_median = median( mode ) ) %>%                              
        ungroup()
      
      dAllTissues_Models_HDI_b_all <- dAllTissues_Models_HDI_b_all %>% 
        bind_rows( ., dAllTissues_Models_HDI_b )
    }
  }
  
  return( dAllTissues_Models_HDI_b_all )
}


plot_figures_3_5_all <- function( m_tisNA = NULL, subgroup_short, 
                                  export = TRUE )
{
  MIN_0 <- get_gr( "MIN_0" )
  f <- get_gr( "f" )
  
  if ( subgroup_short == "entry" ) {
    subgroup <- "Tissue entry\n(rate rel. to blood counts)"
    subgroup_title <- "Tissue entry (rate rel. to blood counts)" }
  if ( subgroup_short == "apoptosis" ) {
    subgroup <- "Apoptosis" 
    subgroup_title <- "Apoptosis" }
  
  d_fig_3_5_rates_all <- get_fig_3_5_rates_all( m_tisNA ) 
  if ( is.null ( d_fig_3_5_rates_all ) ) { return( NULL ) }
    
  d_rates_all <- suppressMessages( 
    d_fig_3_5_rates_all %>% 
    left_join( ., get_pmc()$d_celltypes, by = c( "celltype" = "celltype" ) ) %>% 
    mutate( f.cellstate.group = forcats::fct_recode(
      f.cellstate.group, "Naive/Resting" = "Naive", "CD69+" = "Resident" ) ) %>% 
    group_by( f.celltype, f.tissue.group, f.cellstate.group, 
              par, f.flow.subgroup, f.tissue.group.longer ) %>% 
    summarise( mode = median( mode ) ) %>% 
    ungroup() %>% 
    mutate( mode = ifelse( mode == 0, MIN_0, mode ) )
  )
  
  ymin <- min( d_rates_all$mode * 1000 )
  ymax <- max( d_rates_all$mode * 1000 )

  if ( USE_OLD <- FALSE ) {
  ggfig_4_x <- d_rates_all %>% 
    filter( f.flow.subgroup == subgroup ) %>% 
    ggplot( aes( x = f.celltype, y = mode * 1000 ) ) +
    geom_jitter( aes( shape = f.cellstate.group, colour = f.cellstate.group ),
                 position = position_dodge( 0.3 ), size = 3 * f ) +
    scale_colour_manual( name = "", values = cellstate_color ) +
    scale_y_log10( n.breaks = 8,
                   labels = scales::trans_format( 
                     "log10", scales::math_format( 10^.x ) ),
                   limits = c( ymin, ymax ) ) +
    annotation_logticks( sides = "l" ) +
    labs( title = subgroup_title, 
          y = "Rate (events/1000 cells/day)" ) +     
    facet_grid( . ~ f.tissue.group.longer ) +
    get_gr( "theme_blank" ) +
    theme( legend.title = element_blank(),
           axis.title.y = element_text( hjust = 0.5 ) )
  }
  
  
  ggh4x_scales <- list(
    scale_y_log10( n.breaks = 8,
                   labels = scales::trans_format( 
                     "log10", scales::math_format( 10^.x ) ),
                   limits = c( 10^-1, ymax ) ),
    scale_y_continuous( breaks = c( 0 ),
                        limits = c( -0.01, 0.01 ) ) )
  
  ggfig_4_x_split <- d_rates_all %>% 
    filter( f.flow.subgroup == subgroup ) %>% 
    mutate( is_0 = ifelse( is.na( mode ) | mode == MIN_0, TRUE, FALSE ) ) %>% 
    mutate( mode = ifelse( mode == MIN_0, 0, mode ) ) %>% 
    filter( ! is.na( is_0 ) ) %>% 
    
    ggplot( aes( x = f.celltype, y = mode * 1000 ) ) +
    geom_jitter( aes( shape = f.cellstate.group, colour = f.cellstate.group ),
                 position = position_dodge( 0.3 ), size = 3 * f ) +
    scale_colour_manual( name = "", values = cellstate_color ) +
    scale_shape_manual( name = "", values = cellstate_shape ) +
    labs( title = subgroup_title, y = "Rate (events/1000 cells/day)" ) +     
    facet_grid( is_0 ~ f.tissue.group.longer, 
                scales = "free_y", space = "free_y" ) +
    ggh4x::facetted_pos_scales( y = ggh4x_scales ) +
    ggh4x::force_panelsizes( rows = c( 2, 0.3 ) ) +
    get_gr( "theme_blank" ) +
    theme( legend.title = element_blank(),
           axis.title.y = element_text( hjust = 0.5 ),
           strip.text.y = element_blank() )

  
  
  if ( FALSE ) {
  ggfig_4_x_split2_a <- d_rates_all %>% 
    filter( f.flow.subgroup == subgroup ) %>% 
    bind_rows( ., tibble(
      mode = NA, f.cellstate.group = as.factor( "Naive" ), 
      f.celltype = as.factor( "Tconv" ),
      f.tissue.group.longer = as.factor( "Bone marrow" ) ) ) %>%
    ggplot( aes( x = f.celltype, y = mode * 1000 ) ) +
    geom_jitter( aes( shape = f.cellstate.group, colour = f.cellstate.group ),
                 position = position_dodge( 0.3 ), size = 3 * f ) +
    scale_y_log10( n.breaks = 7,
                   labels = scales::trans_format( 
                     "log10", scales::math_format( 10^.x ) ),
                   limits = c( ymin * 10, ymax ) ) +
  labs( title = subgroup_title, 
        y = "Rate (events/1000 cells/day)" ) +     
    get_gr( "theme_blank" ) +
    theme( legend.title = element_blank(),
           axis.title.y = element_text( hjust = 0.5 ) )
  
  
  ggfig_4_x_split2_b <- d_rates_all %>% 
    filter( f.flow.subgroup == subgroup ) %>% 
    bind_rows( ., tibble(
      mode = NA, f.cellstate.group = as.factor( "Naive" ), 
      f.celltype = as.factor( "Tconv" ),
      f.tissue.group.longer = as.factor( "Bone marrow" ) ) ) %>%
    ggplot( aes( x = f.celltype, y = mode * 1000 ) ) +
    geom_jitter( aes( shape = f.cellstate.group, colour = f.cellstate.group ),
                 position = position_dodge( 0.3 ), size = 3 * f ) +
    scale_y_log10( n.breaks = 3,
                   labels = scales::trans_format( 
                     "log10", scales::math_format( 10^.x ) ),
                   limits = c( ymin, ymin * 10 ) ) +
    labs( title = subgroup_title, 
          y = "Rate (events/1000 cells/day)" ) +     
    get_gr( "theme_blank" ) +
    theme( legend.title = element_blank(),
           axis.title.y = element_text( hjust = 0.5 ) )
  
  ggfig_4_x_split2_ab <- ggpubr::ggarrange( 
    ggfig_4_x_split2_a, ggfig_4_x_split2_b, 
    labels = "",
    ncol = 1, nrow = 2 )  
  }
  
  
  if ( export ) {
    cm_to_in <- get_gr( "cm_to_in" )
    figs_comb_dir <- get_path( "figs_comb", m = m_tisNA )
    if ( ! dir.exists( figs_comb_dir ) ) { dir.create( figs_comb_dir ) }
    fig_fname_start <- get_fig_fname_start( m = m_tisNA, combined = TRUE )
    pdf_name <- sprintf( 
      "%s/%s_Figure_4_%s.pdf", figs_comb_dir, fig_fname_start, subgroup_short )
    ggexport_my( ggfig_4_x_split, filename = pdf_name, 
                 width = 21 * cm_to_in, height = 7 * cm_to_in )
  }
  return( ggfig_4_x_split )
}

plot_PMP_Figure_4 <- function( m_tisNA = m_tisNA )
{
  cm_to_in <- get_gr( "cm_to_in" )
  figs_comb_dir <- get_path( "figs_comb", m = m_tisNA )
  if ( ! dir.exists( figs_comb_dir ) ) { dir.create( figs_comb_dir ) }
  fig_fname_start <- get_fig_fname_start( m = m_tisNA, combined = TRUE )
  pdf_name <- sprintf( "%s/%s_Figure_4AB.pdf", figs_comb_dir, fig_fname_start )
  
  gg_A <- plot_figures_3_5_all( m_tisNA = m_tisNA, subgroup_short = "entry" )
  gg_B <- plot_figures_3_5_all( m_tisNA = m_tisNA, subgroup_short = "apoptosis" )
  gg_AB <- 
    gg_A / gg_B + 
    plot_annotation( tag_levels = "A" ) 
  # https://albert-rapp.de/posts/ggplot2-tips/04_arranging_plots/04_arranging_plots
  
  ggsave( filename = pdf_name, plot = gg_AB,
          width = 21 * cm_to_in, height = 14 * cm_to_in )
  return( gg_AB )
}


if ( FALSE ) {

df <- data.frame(labelx=rep(c('my long label','short'), c(2,26)),
                 labely=rep(c('a','b'), each=14),
                 x=c(letters[1:2],letters[1:26]),
                 y=LETTERS[6:7],
                 i=rnorm(28))
ggplot(df, aes(x,y,color=i)) +
  geom_point() +
  facet_grid(labely~labelx, scales='free_x', space='free_x') +
  force_panelsizes(cols = c(0.3, 2)) +
  theme_bw() 

}



if ( RUN_TEST <- FALSE ) {
  source( "universal_script_setup.r" )
  celltype <- "Tconv"
  source( "testing_model_setup.r" )
  
  plot_PMP_Figure_4( m_tisNA = m_tisNA )
  
  if ( TRUE ) {
  gg_A <- plot_figures_3_5_all( m_tisNA = m_tisNA, subgroup_short = "entry" )
  gg_B <- plot_figures_3_5_all( m_tisNA = m_tisNA, subgroup_short = "apoptosis" )
  gg_AB <- 
    gg_A / gg_B + 
    plot_annotation( tag_levels = "A" ) 
  
  ggsave( filename = "./results/Combined/Figure_4_m0003.pdf", plot = gg_AB,
          width = 21 * 1/2.7, height = 14 * 1/2.7 )
  }
}


