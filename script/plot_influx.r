# plot_influx.r

TESTING_influx <- FALSE

get_prop_sX_orig_s0 <- function( m, n_sim = 10000, sX = "s6", s0 = "s4" ) 
{
  sim_10000 <- sim_combined_dwell( m = m, n_sim = n_sim )
  if ( is.null( sim_10000 ) ) return( NA )
  i_sim_t0_in_sX <-
    sim_10000 %>%
    filter( t == 0, state %in% sX ) %>%
    dplyr::select( i_sim ) %>% unlist() %>% as.integer()
   
  props_orig <-
    sim_10000 %>%
    filter( i_sim %in% i_sim_t0_in_sX ) %>%
    group_by( i_sim ) %>%
    mutate( min_t = min( t ) ) %>%
    ungroup %>%
    filter( t == min_t ) %>%
    summarise( n_s0 = sum( state == s0 ),
               n = n(),
               prop_s0 = n_s0 / n )
  
  return( props_orig$prop_s0 )
}

get_props_activ <- function( m, n_sim = 10000 )
{
  d_orig <- 
    tibble( tissue = get_TISSUES( celltype = m$celltype ), orig_naive = NA )
  for( tissue_j in get_TISSUES( celltype = m$celltype ) ) {
    m$tissue <- tissue_j
    d_orig$orig_naive[ d_orig$tissue == tissue_j ] <- 
      get_prop_sX_orig_s0( m, n_sim = n_sim, sX = "s5", s0 = "s4" )
  }
  return( d_orig )
}


get_props_cd69p <- function( m, n_sim = 10000 )
{
  d_orig <- 
    tibble( tissue = get_TISSUES( celltype = m$celltype ), orig_naive = NA )
  for( tissue_j in get_TISSUES( celltype = m$celltype ) ) {
    m$tissue <- tissue_j
    d_orig$orig_naive[ d_orig$tissue == tissue_j ] <- 
      get_prop_sX_orig_s0( m, n_sim = n_sim, sX = "s6", s0 = "s4" )
  }
  return( d_orig )
}


get_props_activ_cd69p <- function( m, n_sim = 10000 )
{
  d_orig <- 
    tibble( tissue = get_TISSUES( celltype = m$celltype ), orig_naive = NA )
  for( tissue_j in get_TISSUES( celltype = m$celltype ) ) {
    m$tissue <- tissue_j
    d_orig$orig_naive[ d_orig$tissue == tissue_j ] <- 
      get_prop_sX_orig_s0( m, n_sim = n_sim, sX = c( "s5", "s6" ), s0 = "s4" )
  }
  return( d_orig )
}


plot_influx <- function( m, n_sim = 10000 )
{
  m_tisNA <- get_m_tisNA( m = m )
  entry_flow_table_all <- NULL
  for ( tissue in get_pmc()$dTissueAllOrderedGroup$tissue.all.ordered )
  {
    m <- model( m = m_tisNA, tissue = tissue )
    
    if ( dir.exists( get_path( "model_dir", m = m ) ) ) {
      
      entry_flow_table_all_0 <- get_entry_flow_table( 
          ps = m$ps, celltype = m$celltype, tissue = m$tissue,
          model_dir_r = get_path( "model_dir_r", m = m ) ) 
      
      # for incompleted runs - rare issue
      if ( ! is.null( entry_flow_table_all_0 ) ) {
        entry_flow_table_all <-
          entry_flow_table_all_0 %>% 
            mutate( celltype = m$celltype, tissue = tissue ) %>% 
            bind_rows( entry_flow_table_all, . )
      }
    }
  }


  if ( is.null( entry_flow_table_all ) ) { return( NULL ) }

  cm_to_in <- 1 / 2.54   
  f <- cm_to_in * 1.7
  
  if ( TRUE ) {
    THEME <- 
      theme_classic() +
      theme( plot.title = element_text( size = 12 * f, hjust = 0 ), 
             axis.title.x = element_blank(),
             axis.title.y = element_text( size = 12 * f ),  
             
             axis.text.x = element_text( angle = 45, hjust = 1, size = 12 * f ),
             axis.text.y = element_text( size = 12 * f ),  
             
             legend.text = element_text( size = 12 * f ),
             legend.title = element_text( size = 12 * f ),
             
             panel.border = element_blank(),             
             panel.grid.major = element_blank(),         
             panel.grid.minor = element_blank() )
  }
  
  d_tmp <- get_props_activ_cd69p( m, n_sim = n_sim ) %>% 
    left_join( ., get_pmc()$dTissueAllOrderedGroup.f %>% 
                 dplyr::select( tissue.all, f.tissue, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) )
  xintercept1 <- d_tmp$f.tissue[ 
    d_tmp$f.tissue.group == "BoneMarrow" ] %>% unique() %>% length()
  xintercept2 <- d_tmp$f.tissue[ 
    d_tmp$f.tissue.group == "Lymphoid" ] %>% unique() %>% length()
  xintercept3 <- d_tmp$f.tissue[ 
    d_tmp$f.tissue.group == "Non-lymphoid" ] %>% unique() %>% length()
  geom_vline_f.tissues <- geom_vline( xintercept = c( 
      xintercept1 + 0.5, 
      xintercept1 + xintercept2 + 0.5,
      xintercept1 + xintercept2 + xintercept3 + 0.5 ),
      linetype = "dashed" )
  
  
  gg_activ_cd69p <- get_props_activ_cd69p( m, n_sim = n_sim ) %>% 
    left_join( ., get_pmc()$dTissueAllOrderedGroup.f %>% 
                 dplyr::select( tissue.all, f.tissue, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) ) %>% 
    ggplot( aes( x = f.tissue, y = orig_naive ) ) +
    geom_point( size = 3, alpha = 1 ) +
    geom_vline_f.tissues + ylim( 0, 1 ) +
    labs( title = sprintf( 
      "%s, Activated or CD69+ cells entering as naive", m$celltype ),
          y = "Proportion", colour = "Source" ) +     
    THEME
  
  gg_cd69p <- get_props_cd69p( m, n_sim = n_sim ) %>% 
    left_join( ., get_pmc()$dTissueAllOrderedGroup.f %>% 
                 dplyr::select( tissue.all, f.tissue, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) ) %>% 
    ggplot( aes( x = f.tissue, y = orig_naive ) ) +
    geom_point( size = 3, alpha = 1 ) +
    geom_vline_f.tissues + ylim( 0, 1 ) +
    labs( title = sprintf( "%s, CD69+ cells entering as naive", 
                           m$celltype ),
          y = "Proportion", colour = "Source" ) +     
    THEME
  
  gg_activ_check <- get_props_activ( m, n_sim = n_sim ) %>% 
    left_join( ., get_pmc()$dTissueAllOrderedGroup.f %>% 
                 dplyr::select( tissue.all, f.tissue, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) ) %>% 
    ggplot( aes( x = f.tissue, y = orig_naive ) ) +
    geom_point( size = 3, alpha = 1 ) +
    geom_vline_f.tissues + ylim( 0, 1 ) +
    labs( title = sprintf( 
      "%s, Activated cells entering as naive", m$celltype ),
          y = "Proportion", colour = "Source" ) +     
    THEME
  
  dgg_activ <- suppressMessages(
    entry_flow_table_all %>% 
    left_join( ., get_pmc()$dTissueAllOrderedGroup.f %>% 
                 dplyr::select( tissue.all, f.tissue, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) ) %>% 
    filter( population.1 == 5 ) %>% 
    mutate( par_text = ifelse( par == "q25", "Activated from blood", par ) ) %>% 
    mutate( par_text = ifelse( par == "q45", "naive -> activation", par_text ) ) %>%
    mutate( par_text = ifelse( par == "q65", "losing CD69+", par_text ) ) %>% 
    mutate( par_text = as.factor( par_text ) ) %>% 
    mutate( par_group = NA ) %>% 
    mutate( par_group = ifelse( par %in% "q45", "naive", par_group ) ) %>% 
    mutate( par_group = ifelse( par %in% c( "q25", "q65" ), 
                                "activ_cd69p", par_group ) ) %>% 
    group_by( tissue ) %>% 
    mutate( flow_sum = sum( flow ) ) %>% 
    group_by( par_group, f.tissue ) %>% 
    summarise( flow_aggr = sum( flow ) / flow_sum ) %>% 
    ungroup() )
    
  dgg_activ %>% 
    ggplot( aes( x = f.tissue, y = flow_aggr, colour = par_group ) ) +
    geom_point( size = 3, alpha = 1 ) +
    geom_vline_f.tissues + ylim( 0, 1 ) + ylim( 0, 1 ) +
    labs( title = sprintf( "%s, Influx of activated cells", m$celltype ),
          y = "Proportion", colour = "Source" ) +     
    THEME

  if ( PLOT_OLD <- TRUE ) {
    gg_activ_old <- entry_flow_table_all %>% 
      left_join( ., get_pmc()$dTissueAllOrderedGroup.f %>% 
                   dplyr::select( tissue.all, f.tissue, f.tissue.group ), 
                 by = c( "tissue" = "tissue.all" ) ) %>% 
      filter( population.1 == 5 ) %>% 
      mutate( par_text = ifelse( par == "q25", "Activated from blood", par ) ) %>% 
      mutate( par_text = ifelse( par == "q45", "naive -> activation", par_text ) ) %>%
      mutate( par_text = ifelse( par == "q65", "losing CD69+", par_text ) ) %>% 
      mutate( par_text = as.factor( par_text ) ) %>% 
      ggplot( aes( x = f.tissue, y = flow * 1e3, colour = par_text ) ) +
      geom_point( size = 3, alpha = 1 ) + 
      scale_y_log10( n.breaks = 6,
                     limits = c( 0.0001 * 1e3 ,  
                                 1 * 1e3 ),    
                     labels = scales::trans_format( 
                       "log10", scales::math_format( 10^.x ) ) ) +
      geom_vline_f.tissues +
      labs( title = sprintf( "%s, Relative influx of activated cells", m$celltype ),
            y = "Rate (events/1000 cells/day)", colour = "Source" ) +     
      THEME
  }
  
  if ( PLOT_OLD <- TRUE ) {
  gg_cd69p_old <- entry_flow_table_all %>% 
    left_join( ., get_pmc()$dTissueAllOrderedGroup.f %>% 
                 dplyr::select( tissue.all, f.tissue, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) ) %>% 
    filter( population.1 == 6 ) %>% 
    mutate( par_text = ifelse( par == "q36", "CD69+ from blood", par ) ) %>% 
    mutate( par_text = ifelse( par == "q56", "CD69+ differentiation", par_text ) ) %>% 
    mutate( par_text = as.factor( par_text ) ) %>% 
    ggplot( aes( x = f.tissue, y = flow * 1e3, colour = par_text ) ) +
    geom_point( size = 3, alpha = 1 ) +
    scale_y_log10( n.breaks = 6,
                   limits = c( 0.0001 * 1e3, 1 * 1e3 ),    
                   labels = scales::trans_format( 
                     "log10", scales::math_format( 10^.x ) ) ) +
    geom_vline_f.tissues +
    labs( title = sprintf( "%s, Relative influx of CD69+ cells", m$celltype ),
          y = "Rate (events/1000 cells/day)", colour = "Source" ) +     
    THEME
  }
  
  
  fig_filename_start <- get_fig_fname_start( m = m_tisNA )
  figs_dir <- get_path( "figs_m", m = m_tisNA )
  if ( !dir.exists( figs_dir ) ) { dir.create( figs_dir, recursive = TRUE ) }
  
  gg_influx <- ggpubr::ggarrange( 
    gg_cd69p, 
    gg_activ_check,
    labels = c( LETTERS[ 1 : 2 ] ),
    ncol = 2, nrow = 1 )
  pdf_name <- sprintf( "%s/%s_%isims_Figure_5_influx_prop_cd69p.pdf",
                       figs_dir, fig_filename_start, n_sim )
  ggexport_my( gg_influx, filename = pdf_name, 
               width = 10 * 0.8, height = 4.5 * 0.8 )
  
  return( gg_influx )
}

if ( TESTING_influx ) {
  source( "universal_script_setup.r" )
  celltype <- "Treg"
  source( "testing_model_setup.r" )
  
  celltypes <- "Treg"
  celltypes <- c( "Tconv", "CD8", "Treg", "B", "NK" )

  for ( celltype in celltypes ) {
    m_tisNA$celltype <- celltype
    plot_influx( m = m_tisNA, n_sim = 100 )
  }
}

