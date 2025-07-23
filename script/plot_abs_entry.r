# plot_abs_entry.r

plot_abs_entry <- function( celltype = "Tconv", ps_ = ps )
{
  # getwd()
  # celltype <- "Tconv"; ps_ = ps
  total_counts <- read_csv( 
    sprintf( "%s/Total_counts/parabiosis_model_input_%s_counts.csv", 
             ps_$PROCESSED_PATH, celltype ) )

  rel_flows_and_exit <- 
    read_csv( sprintf( "%s/%s/Figures/%s/diag/rel_flows_and_exits.csv", 
                       ps_$RESULTS_PATH, celltype, 
                       paste0( celltype, "_m0003_v01_t1_i1001" ) ) )
  
  total_counts_sum <- sum( total_counts$Median )
  
  d_abs_entry <- 
    rel_flows_and_exit %>% 
    dplyr::select( tissue, al1_Q14, al2_Q25, al3_Q36 ) %>% 
    pivot_longer(
      cols = c(al1_Q14, al2_Q25, al3_Q36),
      names_to = "cellstate",
      values_to = "abs_entry"
    ) %>% 
    mutate( 
      cellstate = ifelse( cellstate == "al1_Q14", "Naive/Resting", cellstate ),
      cellstate = ifelse( cellstate == "al2_Q25", "Antigen-experienced", cellstate ),
      cellstate = ifelse( cellstate == "al3_Q36", "CD69+", cellstate ) ) %>% 
    mutate( cellstate = factor( 
      cellstate, c( "Naive/Resting", "Antigen-experienced", "CD69+" ) ) ) %>% 
    left_join( ., get_pmc()$dTissueAllOrderedGroup.f, 
               by = c( "tissue" = "tissue.all" ) )
  
  xintercept1 <- d_abs_entry$f.tissue[ 
    d_abs_entry$f.tissue.group == "BoneMarrow" ] %>% unique() %>% length()
  xintercept2 <- d_abs_entry$f.tissue[ 
    d_abs_entry$f.tissue.group == "Lymphoid" ] %>% unique() %>% length()
  xintercept3 <- d_abs_entry$f.tissue[ 
    d_abs_entry$f.tissue.group == "Non-lymphoid" ] %>% unique() %>% length()
  
  gg_abs_entry <- 
    d_abs_entry %>% 
    ggplot( aes( x = f.tissue, y = abs_entry * total_counts_sum, 
                 colour = cellstate, shape = cellstate ) ) +
    geom_jitter( position = position_dodge( 0.3 ) ) +
    geom_vline( xintercept = 
                  c( xintercept1 + 0.5, 
                     xintercept1 + xintercept2 + 0.5,
                     xintercept1 + xintercept2 + xintercept3 + 0.5 ) ) +
    scale_colour_manual( name = "", values = cellstate_color ) +
    scale_shape_manual( name = "", values = cellstate_shape ) +
    scale_x_discrete( labels = c( "BoneMarrow" = "Bone marrow" ) ) +
    labs( title = celltype, y = "Cumulative entry (cells/day)" ) +
    scale_y_log10( labels = scales::trans_format( 
      "log10", scales::math_format( 10^.x ) ) ) + 
    annotation_logticks( sides = "l" ) +
    get_gr( "theme_blank" ) + 
    theme( axis.title.y = element_text( hjust = 0.5 ) ); gg_abs_entry
  
  return( gg_abs_entry )
}


if ( RUN_TEST <- FALSE ) {
  library( patchwork )
  FILL_IN_MY_PROJECT_PATH <- ""
  MY_PROJECT_PATH <- ifelse( nchar( FILL_IN_MY_PROJECT_PATH ) != 0,
                             FILL_IN_PROJECT_PATH, getwd() )
  source( sprintf( "%s/SET_PATHS.r", MY_PROJECT_PATH ) )
  PS <- SET_PATHS( PROJECT_PATH = MY_PROJECT_PATH )
  ps <- PS;

  celltypes <- c( "CD8", "Tconv", "Treg" )
  l_gg_en <- list() # ; N_sim <- 10000; X_max_cut <- 1500
  for ( celltype in celltypes ) {
    i <- which( celltypes == celltype )
    # for ( state in states ) {
      # d_probs <- get_d_probs( celltype = celltype, state = state )
      # d_sim <- sim_lifespan( n_sim = N_sim, d_ = d_probs )
    l_gg_en[[ i ]] <- plot_abs_entry( celltype = celltype )
      # plot_lifespan( celltype = celltype, state = state, d_sim_ = d_sim, 
      #                x_max_cut = X_max_cut )
    # }
  }
  gg_en_3 <- l_gg_en[[ 1 ]] / l_gg_en[[ 2 ]] / l_gg_en[[ 3 ]] + 
    plot_annotation( tag_levels = "A" ); gg_en_3
  ggsave( filename = sprintf( 
    "./results/Combined/abs_entry_%s.pdf",
    "CD8-Tconv-Treg" ),
    plot = gg_en_3, width = 21 / 2.7, height = 29.7 / 2.7 )
}
