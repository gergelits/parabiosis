# sim_tissue_state_model.r

sim_combined_dwell <- function( m, n_sim = 1000, incl_past = TRUE,
                                Q_from_HDI = TRUE )
{
  chain <- get_max_chain( m )
  source( sprintf( "%s/%s", m$ps$CODE_PATH, "compare_qij_and_mode.r"), 
          local = TRUE ) 
  chain_dir <- get_path( "chain_dir", m = m, sel_type = "max",
                         sel_chains = get_max_chain( m ) )
  fig_fname_start <- get_fig_fname_start( m = m )
  figs_dwell_dir <- get_path( "figs_m_dwell", m = m )
  sim_file <- sprintf( "%s/%s_%s_%isims_steps_%s.csv", 
                       figs_dwell_dir, fig_fname_start, 
                       tolower( m$tissue ), n_sim,
                       ifelse( incl_past, "past_fut", "fut" ) )
  sim_done <- file.exists( sim_file )
  
  # chain <- get_max_chain( m )    
  parabio_fit_HDI <- get_parabio_fit_HDI( m = m, chain = chain )
  if ( Q_from_HDI & ( ! is.null( parabio_fit_HDI ) ) ) { 
    
    Q <- HDI_table_to_Q_matrix( parabio_fit_HDI ) } else
      Q <- read_Q_matrix( chain_dir )
  
  if ( is.null( Q ) ) { return( NULL ) }
  
  if ( ! sim_done ) {
    al_ <- get_al_( m = m )
    source( sprintf( "%s/%s", m$ps$CODE_PATH, "tissue_states_model.r" ), 
            local = TRUE )
    m456_ex <- get_tissue_states_model_exit( Q = Q, al_ = al_ )
    m456_en <- get_tissue_states_model_entry( Q = Q, al_ = al_ )
    sim_f <- sim_tissue_states_dwell_time_exit( n_sim = n_sim, m456 = m456_ex )
    
    if ( incl_past ) {
      sim <- sim_tissue_states_past_future( m456_en = m456_en, 
                                            sim_future = sim_f )
    } else { sim <- sim_f }
    sim %>% write_csv( sim_file )
  } else {
    sim <- read_csv( sim_file )
  }
  
  return( sim )
}


plot_combined_dwell <- function( m, sim, ylim_max = NA, export = TRUE,
                                 legend = "out", legend_show = TRUE,
                                 log_time = FALSE )
{
  if ( legend == "out" ) { theme_legend <- theme( legend.position = "right" ) }
  if ( legend == "inplot" ) { 
    theme_legend <- theme( legend.position = c( 0.8, 0.8 ) ) }
  if ( !legend_show ) { theme_legend <- theme( legend.position = "none" ) }
  
  sim_steps_ordered <- sim %>%
    group_by( i_sim ) %>% 
    summarise( total_t = sum( time ) ) %>% 
    mutate( i_sim_order = rank( -total_t ) )
  
  # output:
  # paths
  fig_fname_start <- get_fig_fname_start( m = m )
  figs_dwell_dir <- get_path( "figs_m_dwell", m = m )
  
  # table / text
  prop_exit <- sim %>% 
    filter( occur == "exit" ) %>% 
    group_by( state ) %>% 
    summarise( n = n() ) %>% 
    mutate( freq = n / sum( n ) )
  
  # Check to have all "state" in table even if it equals 0.
  prop_exit <- 
    prop_exit %>% 
    bind_rows( tibble( 
      state = c( "blood_s5", "blood_s6", "die_s5", "die_s6", "exit_s4" ),
      n = rep( 0, 5 ),
      freq = rep( 0, 5 ) ) ) %>% 
    group_by( state ) %>% 
    filter( n == max( n ) ) %>% unique() %>% 
    ungroup() %>% arrange( state )
  
  if ( export ) {
    n_sim <- max( sim$i_sim )
    tibble( tissue = tolower( m$tissue ),
            median = median( sim_steps_ordered$total_t ) ) %>% 
      write_csv( sprintf( "%s/%s_%s_%isims_median_dwell.csv", figs_dwell_dir, 
                          fig_fname_start, tolower( m$tissue ), n_sim ) )
    
    prop_exit %>% 
      write_csv( sprintf( "%s/%s_%s_%isims_prop_exit.csv", figs_dwell_dir, 
                          fig_fname_start, tolower( m$tissue ), n_sim ) )
  }
    
  # figures
  cm_to_in <- 1 / 2.54   
  f <- cm_to_in * 1.7
  pal <- hue_pal()( 3 )
  pal_vector <- pal[ c( rep( c( 2, 3 ), 100 ), 1 ) ]
  names( pal_vector ) <- sprintf( 
    "o%s", stringr::str_pad( c( 1 : ( length( pal_vector ) - 1 ), 0 ), 
                             3, pad = "0" ) )

  sim_steps_ordered <- 
    sim_steps_ordered %>% filter( !is.na( total_t ) )
  
  q50 <- round( quantile( sim_steps_ordered$total_t, 0.50 ) )
  q95 <- round( quantile( sim_steps_ordered$total_t, 0.95 ) )
  
  # Plot only subset of 1000 as readable. While estimated from more.
  sim_sub <- sim %>% 
    filter( i_sim <= 1000 )
  sim_steps_ordered_sub <- sim_sub %>%
    group_by( i_sim ) %>% 
    summarise( total_t = sum( time ) ) %>% 
    mutate( i_sim_order = rank( -total_t ) )
  n_sim <- max( sim_sub$i_sim )
  
  sim_sub_to_plot <- sim_sub %>%
    left_join( ., sim_steps_ordered_sub %>% dplyr::select( i_sim, i_sim_order ),
               by = c( "i_sim" = "i_sim" ) ) %>%
    mutate( occur = factor( 
      occur, levels = sort( unique( occur ), decreasing = TRUE ) ) ) %>% 
    group_by( i_sim ) %>% 
    mutate( cum_time = cumsum( time ) ) %>% 
    group_by( i_sim, occur ) %>% 
    mutate( time_cut = ifelse( 
      ( ! is.na( ylim_max ) & ( occur != "exit" ) & ( cum_time > ylim_max ) ), 
      max( ylim_max - ( cum_time - time ), 0 ), time ) ) %>% 
    mutate( time = time_cut )
  

  q95_shft <- ifelse( ( q95 < ylim_max & q95 > ylim_max * 0.75 ), 
                      q95 - ylim_max * 0.21, q95 )
  
    ggfig_steps <- sim_sub_to_plot %>% 
    ggplot( aes( x = i_sim_order, y = time, fill = occur ) ) +
    geom_bar( position = "stack", stat = "identity" ) + 
    labs( title = sprintf( "%s - %s", m$celltype, m$tissue ), 
          x = "Sampled cell index", 
          y = "Dwell time (days)",
          fill = "Cell state" ) + 
    ylim( 0, ylim_max ) + 
    theme_legend +
    coord_flip() +
    
    geom_hline( yintercept = q50, linetype = "dashed" ) + 
    annotate( x = 0.85 * n_sim, y = min( ylim_max * 0.75, q50, na.rm = TRUE ), 
              hjust = -0.1, geom = "text", colour = "black", 
              label = sprintf( "q50=%s", q50 %>% as.character() ) ) +
    
    geom_hline( yintercept = q95, linetype = "dashed", colour = "darkred" ) + 
    
    annotate( x = 0.15 * n_sim, 
              y = min( ylim_max * 0.75, q95_shft, na.rm = TRUE ), 
              hjust = -0.1, geom = "text", colour = "darkred", 
              label = sprintf( "q95=%s", q95 %>% as.character() ) ) + 
    scale_fill_manual( values = pal_vector,
                       guide = guide_legend( reverse = TRUE ),
                       breaks = c( "o000", "o001", "o002" ),
                       labels = c( "naive", "activated", "CD69+" ) )
  
  if ( export ) {
    pdf_name <- sprintf( "%s/%s_%s_%isims_dwell_steps.pdf", figs_dwell_dir, 
                         fig_fname_start, tolower( m$tissue ), n_sim )
    ggexport_my( ggfig_steps, filename = pdf_name,
                 width = 21 * cm_to_in, height = 6 * cm_to_in )
  }

  return( ggfig_steps )
}


get_d_probs <- function( celltype, state, ps_ = ps )
{
  model_id<-sprintf("%s_m0003_v01_t1_i1001",celltype)
  f_rel_flows <- file.path( ps_$RESULTS_PATH, celltype, "Figures", model_id, 
                            "diag", "rel_flows_and_exits.csv" )
  
  f_dwell <- file.path( ps_$RESULTS_PATH, celltype, "Figures", model_id, 
                        sprintf( "%s_dwell_summary.csv", model_id ) )
  
  d_dwell <- read_csv( f_dwell ) %>% 
    dplyr::select( tissue, activ, cd69p ) %>% 
    filter( ! ( is.na( tissue ) ) & 
              tissue != "Median" ); d_dwell
  
  if ( state == "cd69p" ) {
    d_probs <- read_csv( f_rel_flows ) %>% 
      dplyr::select( tissue, al3_Q36, blood_s6, die_s6 ) %>% 
      mutate( al3_Q36_sum = sum( al3_Q36 ) ) %>% 
      mutate( tissue_i = al3_Q36 / al3_Q36_sum ) %>% 
      mutate( P6_sum = blood_s6 + die_s6,
              die_i = die_s6 / P6_sum, 
              blood_i = blood_s6 / P6_sum ) %>% 
      left_join( ., d_dwell, by = "tissue" ) %>% 
      dplyr::rename( dwell = cd69p ); d_probs
  } else {
    d_probs <- read_csv( f_rel_flows ) %>% 
      dplyr::select( tissue, al2_Q25, blood_s5, die_s5 ) %>% 
      mutate( al2_Q25_sum = sum( al2_Q25 ) ) %>% 
      mutate( tissue_i = al2_Q25 / al2_Q25_sum ) %>% 
      mutate( P5_sum = blood_s5 + die_s5,
              die_i = die_s5 / P5_sum, 
              blood_i = blood_s5 / P5_sum ) %>% 
      left_join( ., d_dwell, by = "tissue" ) %>% 
      dplyr::rename( dwell = activ ); d_probs
  }

  return( d_probs )
}

sim_lifespan <- function( m = NULL, n_sim = 1000, 
                          d_ = d_probs, s_init = "blood" )
{
  sim_steps <- NULL
  for ( i in 1 : n_sim ) {
    if ( i %% 200 == 0 ) { cat( "sim ", i, "\n" ) }
    t <- 1
    if ( s_init == "blood" ) {
      s_curr <- sample( x = d_$tissue, size = 1, prob = d_$tissue_i )
    } else {
      s_curr <- s_init
    }
    cell_next <- "alive"
    while ( cell_next != "die" ) {
      curr_dwell_time <- 
        rexp( n = 1, rate = 1 / d_$dwell[ d_$tissue == s_curr ] )
    
      
      sim_step_curr <- tibble( 
        i_sim = i, 
        tissue = s_curr,
        t = t,
        time = curr_dwell_time )
      sim_steps <- bind_rows( sim_steps, sim_step_curr )
      
      cell_next <- sample( x = c( "in", "die", "blood" ), size = 1,
                             prob = c( 0,
                                       d_$die_i[ d_$tissue == s_curr ],
                                       d_$blood_i[ d_$tissue == s_curr ] ) )
      if ( cell_next != "die" ) {
        s_curr <- sample( x = d_$tissue, size = 1, prob = d_$tissue_i )
        t <- t + 1
      }
    }
  }
  return( sim_steps )
}


plot_lifespan <- function( celltype = "Treg", state = "cd69p", d_sim_ = d_sim,
                           x_max_cut = 1000 )
{
  binw <- 50
  d_sim_summ_ <- d_sim %>% 
    group_by( i_sim ) %>% 
    summarise( time_sum = sum( time ) ) %>% 
    mutate( time_sum = ifelse( time_sum >= x_max_cut, x_max_cut, time_sum ) )
  
  
  time_sum_med <- d_sim_summ_$time_sum %>% median()
  plot_build <- ggplot_build(
    ggplot( d_sim_summ_, aes( x = time_sum ) ) + 
      geom_histogram( binwidth = binw ) )
  
  max_count <- plot_build$data[[ 1 ]]$count %>% max()
  bin_breaks <- plot_build$data[[ 1 ]]$x
  bin_breaks <- seq( 0, x_max_cut, 100 )
  bin_breaks_lbl <- ifelse( bin_breaks >= x_max_cut, 
                            sprintf( "%i+", x_max_cut ), bin_breaks )  
    
  gg <- d_sim_summ_ %>% 
    ggplot( aes( x = time_sum ) ) +
    labs( x = "Cell lifespan (days)", y = "# simulations",
          ) +
    geom_histogram( binwidth = binw ) +
    scale_x_continuous( breaks = bin_breaks, labels = bin_breaks_lbl ) +
    geom_vline( xintercept = time_sum_med, color = "red" ) + 
    annotate( "text", color = "red",
              x = time_sum_med * 1.1, y = max_count * 0.85, hjust = 0,
              label = paste0( "q50 = ", round( time_sum_med ), " days" ) )
  print( gg )
  
  ggsave( filename = sprintf( "./results/Combined/lifespan_%s_%s_%ssim_max%s.pdf",
                              celltype, state, 
                              d_sim_summ_$i_sim %>% max %>% as.character,
                              x_max_cut ),
          plot = gg, width = 10, height = 6 )
  return( gg )
}


if ( RUN_TEST <- FALSE ) {
  
  FILL_IN_MY_PROJECT_PATH <- ""
  MY_PROJECT_PATH <- ifelse( nchar( FILL_IN_MY_PROJECT_PATH ) != 0,
                             FILL_IN_PROJECT_PATH, getwd() )
  source( sprintf( "%s/SET_PATHS.r", MY_PROJECT_PATH ) )
  library( tidyverse )
  library( patchwork )
  PS <- SET_PATHS( PROJECT_PATH = MY_PROJECT_PATH )

  ps <- PS;
  celltypes <- c( "CD8", "Tconv", "Treg" )
  states <- c( "cd69p" )
  # states <- c( "activ", "cd69p" )
  
  l_gg <- list(); N_sim <- 10000; X_max_cut <- 1000
  for ( celltype in celltypes ) {
    i <- which( celltypes == celltype )
    for ( state in states ) {
      d_probs <- get_d_probs( celltype = celltype, state = state )
      d_sim <- sim_lifespan( n_sim = N_sim, d_ = d_probs )
      l_gg[[ i ]] <- 
        plot_lifespan( celltype = celltype, state = state, d_sim_ = d_sim, 
                       x_max_cut = X_max_cut )
      plot_lifespan( celltype = celltype, state = state, d_sim_ = d_sim, 
                     x_max_cut = X_max_cut )
    }
  }
  gg_3 <- l_gg[[ 1 ]] / l_gg[[ 2 ]] / l_gg[[ 3 ]] + 
    plot_annotation( tag_levels = "A" ); gg_3
  ggsave( filename = sprintf( 
    "./results/Combined/lifespan_%s_%s_%isim_max%s.pdf",
    "CD8-Tconv-Treg", state, N_sim, X_max_cut ),
          plot = gg_3, width = 21 / 2, height = 29.7 / 2 )
}

