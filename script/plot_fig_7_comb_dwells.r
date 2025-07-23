# plot_fig_7_comb_dwells.r

plot_extra_comb_dwell <- function( l_gg_all, m, celltype, tissue, 
                                   n_sim = 10000, ylim_max = 500 )
{
  max_k <- length( l_gg_all )
  m <- model( m = m, celltype = celltype, tissue = tissue )
  sim_10000 <- sim_combined_dwell( m = m, n_sim = n_sim )
  if ( ! is.null( sim_10000 ) ) { 
    l_gg_all[[ max_k + 1 ]] <- 
      plot_combined_dwell( m = m, sim = sim_10000, ylim_max = ylim_max, 
                           export = FALSE, legend_show = FALSE ) 
  } else { 
    l_gg_all[[ max_k + 1 ]] <- ggplot() + theme_void()
    }
  return( l_gg_all )
}

plot_fig_7_comb_dwells <- function( m_tisNA, n_sim = 10000, export = TRUE,
                                    ylim_max = 50 )
{
  l_gg_all <- list()
  celltypes <- c( "CD8", "Tconv", "Treg", "B", "NK" )
  tissues <- c( "Spleen", "Brain", "Liver", "Lung", "WAT", "IEL" )
  model_name <- m_tisNA$model_name
  
  m_tisNA$model_id <- NULL  
  
  stan_pars_v_0003 <- "v01"
  m_tisNA_0003 <- model( 
    model = m_tisNA, 
    model_name = model_name,
    stan_pars_v = stan_pars_v_0003, 
    model_ver = sprintf( "%s_t%s", stan_pars_v_0003, N_TISSUES = 1 ) )

  for ( celltype in celltypes ) {
    
    m_tisNA <- m_tisNA_0003
    m_tisNA$celltype <- celltype
    
    j <- which( celltype == celltypes )
    for ( tissue in tissues ) {
      i <- which( tissue == tissues )
      k <- ( j - 1 ) * length( tissues ) + i
      l_gg_all[[ k ]] <- ggplot() + theme_void()  # patchwork::plot_spacer()
      m <- model( m = m_tisNA, tissue = tissue )
      sim_10000 <- sim_combined_dwell( m = m, n_sim = n_sim )
      if ( ! is.null( sim_10000 ) ) {
        l_gg_all[[ k ]] <- plot_combined_dwell( 
          m = m, sim = sim_10000, ylim_max = ylim_max, export = FALSE,
          legend = "inplot", legend_show = ( k == 1 ) )
          # legend = "inplot", legend_show = ( k == 6 ) )
      } 
    }
  }
  
  
  d_extra_row <- tibble( 
    celltype = c( "CD8", "Tconv", "Tconv", "Treg", "NK", "NK" ),
    tissue   = c( "IEL", "Liver", "IEL", "IEL",  "Liver",  "IEL" ) )
  for ( i in 1 : nrow( d_extra_row ) ) {
    l_gg_all <- plot_extra_comb_dwell( 
      l_gg_all = l_gg_all, m = m_tisNA, 
      celltype = d_extra_row$celltype[ i ], tissue = d_extra_row$tissue[ i ],
      n_sim = n_sim )  
  }
  
  ncol = length( tissues )
  nrow = ceiling( length( l_gg_all ) / length( tissues ) )
  labels = matrix( "", ncol = ncol, nrow = nrow )
  labels[ , 1 ] <- LETTERS[ 1 : nrow ]
  labels <- labels %>% t() %>% unlist() %>% as.vector()
  gg_examples_comb_dwells <- ggpubr::ggarrange( 
    plotlist = l_gg_all, ncol = ncol, nrow = nrow, labels = labels )
  
  if ( export ) {
    cm_to_in <- get_gr( "cm_to_in" )
    figs_comb_dir <- get_path( "figs_comb", m = m_tisNA )
    if ( ! dir.exists( figs_comb_dir ) ) { dir.create( figs_comb_dir ) }
    fig_fname_start <- get_fig_fname_start( m = m_tisNA, combined = TRUE ) 
    pdf_name <- sprintf( 
      "%s/%s_%isims_Figure_7_comb_dwell_examples_%s.pdf", 
      figs_comb_dir, fig_fname_start, n_sim, 
      ifelse( !is.na( ylim_max ), as.character( ylim_max ), "rel" ) )
    ggexport_my( gg_examples_comb_dwells, filename = pdf_name, 
                 width = 29.7 * cm_to_in * 2, height = 21 * cm_to_in * 2 )
  }
  
  return( gg_examples_comb_dwells )
}

# Need to be run from outside to avoid deep recursion:
if ( RUN_TEST <- FALSE ) {
  source( "universal_script_setup.r" )
  celltype <- "Tconv"
  source( "testing_model_setup.r" )
  gg <- plot_fig_7_comb_dwells( m_tisNA = m_tisNA, n_sim = 1000, ylim_max = 50 )
}
