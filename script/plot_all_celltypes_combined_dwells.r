# plot_all_celltypes_combined_dwells.r

get_d_res <- function( m_tisNA )
{
  celltypes <- c( "Tconv", "CD8", "Treg", "B", "NK" ) 
  
  d_res <- tibble(
    celltype = celltypes,
    mcmc     = rep( m_tisNA$mcmc_pars$id, length( celltypes ) ),
    model_id = rep( m_tisNA$model_id,     length( celltypes ) )
  )
  
  return( d_res )
}

read_medians <- function( d_res, m_tisNA, n_sim = 10000 )
{
  d_medians <- tibble( tissue = NULL, median = NULL )
  for ( i in 1 : nrow( d_res ) ) {
    celltype <- d_res$celltype[ i ]
    mcmc <- d_res$mcmc[ i ]
    model_id <- d_res$model_id[ i ]
    
    for ( tissue_j in get_TISSUES( celltype = celltype ) )
    {
      stan_pars_v <- m_tisNA$stan_pars_v; ps <- m_tisNA$ps
      m <- model( ps = ps,
                  celltype = celltype,
                  tissue = tissue_j, 
                  mcmc_pars = get_mcmc_pars( mcmc_pars_v = mcmc, ps = ps ),
                  stan_pars_v = stan_pars_v,
                  model_name = sprintf( "model_%s", model_id ),
                  model_ver = sprintf( "%s_t%s", stan_pars_v, N_TISSUES = 1 ),
                  max_lp_chain = 1 )
      
      fig_fname_start <- get_fig_fname_start( m = m )
      figs_dwell_dir <- get_path( "figs_m_dwell", m = m )
      
      i_j_file <- sprintf( "%s/%s_%s_%ssims_median_dwell.csv", figs_dwell_dir, 
        fig_fname_start, tolower( m$tissue ), n_sim ) 
      if ( file.exists( i_j_file ) ) {
        d_medians <- read_csv( i_j_file ) %>% 
          mutate( celltype = celltype ) %>% 
          bind_rows( d_medians, . )
          
      }
      # else: nothing is added to d_dwell_cd69p for this celltype
    }
  }
  return( d_medians )
}

read_dwell_cd69p <- function( d_res, m_tisNA )
{
  d_dwell_cd69p <- NULL
  for ( i in 1 : nrow( d_res ) ) {
    celltype <- d_res$celltype[ i ]
    mcmc <- d_res$mcmc[ i ]
    model_id <- d_res$model_id[ i ]
    
    stan_pars_v <- m_tisNA$stan_pars_v; ps <- m_tisNA$ps
    m_tisNA <- model( ps = ps,
                celltype = celltype,
                mcmc_pars = get_mcmc_pars( mcmc_pars_v = mcmc, ps = ps ),
                stan_pars_v = stan_pars_v,
                model_name = sprintf( "model_%s", model_id ),
                model_ver = sprintf( "%s_t%s", stan_pars_v, N_TISSUES = 1 ),
                max_lp_chain = 1 )
    
    fig_fname_start <- get_fig_fname_start( m = m_tisNA )
    figs_dir <- get_path( "figs_m", m = m_tisNA )
    
    csv_name <- sprintf( "%s/%s_dwell_CrI.csv", figs_dir, fig_fname_start )
    
    if ( file.exists( csv_name ) ) {
      csv_tmp <- read_load_multi( fun = "read_csv", full_path = csv_name )
    
      d_dwell_cd69p <- 
        csv_tmp %>% 
        filter( cellstate == "CD69+" ) %>% 
        dplyr::select( celltype, tissue, dwell_time_mode_DAYS,	
                       dwell_time_hdi_80_upper,	dwell_time_hdi_80_lower) %>% 
        bind_rows( d_dwell_cd69p, . )
    }
    # else: nothing is added to d_dwell_cd69p for this celltype
  } 

  return( d_dwell_cd69p )
}


plot_all_dwell <- function( d_dwells = d_dwells_combined, m_tisNA, 
                            n_sim = 10000, dwell_max = 120 )
{
  d_dwells_groups <- d_dwells %>% 
    mutate( tissue = tolower( tissue ) ) %>% 
    left_join( ., get_pmc()$dTissueAllOrderedGroup.f %>% 
                 dplyr::select( tissue.all, f.tissue, f.tissue.group ) %>% 
                 mutate( lower.tissue.all = tolower( tissue.all ) ), 
               by = c( "tissue" = "lower.tissue.all" ) ) %>% 
    
    left_join( ., get_pmc()$d_celltypes, by = c( "celltype" = "celltype" ) ) %>% 
    mutate( celltype = f.celltype ) %>% 
    filter( f.tissue.group != "BoneMarrow" ) %>% 
    mutate( text_outlier = ifelse( 
      dwell > dwell_max - 20, sprintf( "%i d", dwell %>% round() ), "" ) )
    
  
  xintercept1 <- d_dwells_groups$f.tissue[ 
    d_dwells_groups$f.tissue.group == "BoneMarrow" ] %>% unique() %>% length()
  xintercept2 <- d_dwells_groups$f.tissue[ 
    d_dwells_groups$f.tissue.group == "Lymphoid" ] %>% unique() %>% length()
  xintercept3 <- d_dwells_groups$f.tissue[ 
    d_dwells_groups$f.tissue.group == "Non-lymphoid" ] %>% unique() %>% length()
  
  cm_to_in <- 1 / 2.54   
  f <- cm_to_in * 1.7
  
  THEME <- theme_classic() +
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
  
  gg_all_dwell <- d_dwells_groups %>%
    mutate( dwell = ifelse( dwell > dwell_max, dwell_max, dwell ) ) %>%
    ggplot( aes( x = f.tissue, y = dwell, fill = celltype ) ) +
    geom_bar( stat = "identity", position = "dodge" ) +
    geom_text( aes( label = text_outlier, x = f.tissue, y = dwell_max + 2 ),
               hjust = 0.5, vjust = 0, position = position_dodge( .9 ), size = 3 ) +
    geom_rect( aes( xmin = f.tissue %>% as.numeric() - 1.45, 
                    xmax = f.tissue %>% as.numeric() - .55, 
                    ymin = dwell_max - 20, ymax = dwell_max - 10 ),
               size = 0.3,
               fill = "white" ) +
    geom_vline( xintercept = c( xintercept1 + 0.5, 
                                xintercept1 + xintercept2 + 0.5,
                                xintercept1 + xintercept2 + xintercept3 + 0.5 ),
                linetype = "dashed" ) +
    labs( y = "Median Dwell Time (days)", fill = "Cell type" ) +    
    scale_y_continuous( breaks = seq( from = 0, to = dwell_max - 20, by = 20 ),
                        limits = c( 0, dwell_max + 12 ) ) +   # , labels = y_labels)
    THEME; gg_all_dwell
  if ( FALSE ) {
  # NOTE:
  # Or use split y axis as in 
  # https://rstudio-pubs-static.s3.amazonaws.com/513019_8f16079f34f04ebf98d86e8dcbcc9c21.html
  # Not working well for barcharts and my situation
  gg_all_dwell_split <- d_dwells_groups %>%
    mutate( facet_bin = dwell > dwell_max ) %>% 
    
    ggplot( aes( x = f.tissue, y = dwell, fill = celltype ) ) +
    geom_bar( stat = "identity", position = "dodge" ) +
    geom_text( aes( label = text_outlier, x = f.tissue, y = dwell ),
               vjust = 1.5, position = position_dodge( .9 ), size = 3 ) +
    geom_vline( xintercept = c( xintercept1 + 0.5, 
                                xintercept1 + xintercept2 + 0.5,
                                xintercept1 + xintercept2 + xintercept3 + 0.5 ),
                linetype = "dashed" ) +
    labs( y = "Median Dwell Time (days)", fill = "Cell type" ) +    
    facet_grid( facet_bin ~ ., scale = "free_y" ) +
    THEME
  }
  
  figs_comb_dir <- get_path( "figs_comb", m = m_tisNA )
  if ( ! dir.exists( figs_comb_dir ) ) { dir.create( figs_comb_dir ) }
  fig_fname_start <- get_fig_fname_start( m = m_tisNA, combined = TRUE ) 
  
  if ( unique( d_dwells$type ) == "combined" ) { 
    pdf_name <- sprintf( "%s/%s_%isims_Figure_6C_comb_dwell_medians.pdf", 
                         figs_comb_dir, fig_fname_start, n_sim )
  }
  if ( unique( d_dwells$type ) == "cd69p" ) { 
    pdf_name <- sprintf( "%s/%s_Figure_6D_cd69p_dwell_medians.pdf", 
                         figs_comb_dir, fig_fname_start )
  }
  
  ggexport_my( gg_all_dwell, filename = pdf_name, 
               width = 29.7 * cm_to_in, height = 21 / 3 * cm_to_in )
  
  return( gg_all_dwell )
}


plot_all_celltypes_dwells <- function( m_tisNA, n_sim )
{
  d_res <- get_d_res( m_tisNA = m_tisNA )
  
  d_dwells_cd69p <- read_dwell_cd69p( d_res = d_res, m_tisNA = m_tisNA ) %>% 
    dplyr::rename( dwell = dwell_time_mode_DAYS ) %>% 
    mutate( type = "cd69p" )
  
  d_dwells_combined <- read_medians( d_res = d_res, m_tisNA = m_tisNA, 
                                     n_sim = n_sim ) %>% 
    dplyr::rename( dwell = median ) %>% 
    mutate( type = "combined" )
  
  gg_combined <- plot_all_dwell( d_dwells = d_dwells_combined, 
                                 m_tisNA = m_tisNA, n_sim = n_sim )
  
  gg_cd69p <- plot_all_dwell( d_dwells = d_dwells_cd69p, m_tisNA = m_tisNA )
}


if ( FALSE ) { source( "universal_script_setup.r" ) }

if ( FALSE ) { 
  source( "universal_script_setup.r" )
  celltype <- "NK"
  source( "testing_model_setup.r" )
  plot_all_celltypes_dwells( m_tisNA = m_tisNA, n_sim = 100 )
}
