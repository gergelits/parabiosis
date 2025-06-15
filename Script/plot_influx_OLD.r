# plot_influx.r
if ( FALSE ) {
if ( FALSE ) {
source( "universal_script_setup.r" )
celltype <- "CD8"
source( "testing_model_setup.r" )

celltypes <- c( "Tconv", "CD8", "Treg", "B", "NK" )
celltypes <- c( "NK" ) 

plot_influx <- function( m )
{
  m_tisNA <- get_m_tisNA( m = m )
  entry_flow_table_all <- NULL
  for ( tissue in get_pmc()$dTissueAllOrderedGroup$tissue.all.ordered )
  {
    m <- model( m = m_tisNA, tissue = tissue )
    
    if ( dir.exists( get_path( "model_dir", m = m ) ) ) {
      entry_flow_table_all <- get_entry_flow_table( 
          ps = m$ps, celltype = m$celltype, tissue = m$tissue,
          model_dir_r = get_path( "model_dir_r", m = m ) ) %>% 
        mutate( celltype = m$celltype, tissue = tissue ) %>% 
        bind_rows( entry_flow_table_all, . )
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
           axis.title.y = element_text( size = 12 * f, hjust = 1 ),  
           
           axis.text.x = element_text( angle = 45, hjust = 1, size = 12 * f ),
           axis.text.y = element_text( size = 12 * f ),  
           
           legend.text = element_text( size = 12 * f ),
           legend.title = element_text( size = 12 * f ),
           
           panel.border = element_blank(),             
           panel.grid.major = element_blank(),         
           panel.grid.minor = element_blank() )
  }
  
  gg_activ <- entry_flow_table_all %>% 
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
    
    labs( title = sprintf( "%s, Relative influx of acitvated cells", m$celltype ),
          y = "Rate (events/1000 cells/day)", colour = "Source" ) +     
    THEME
  
  
  
  gg_cd69p <- entry_flow_table_all %>% 
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
                 limits = c( 0.0001 * 1e3 ,  
                             1 * 1e3 ),    
                 labels = scales::trans_format( 
                   "log10", scales::math_format( 10^.x ) ) ) +
    
    labs( title = sprintf( "%s, Relative influx of CD69+ cells", m$celltype ),
          y = "Rate (events/1000 cells/day)", colour = "Source" ) +     
    THEME
  
  fig_filename_start <- get_fig_fname_start( m = m_tisNA )
  
  figs_dir <- get_path( "figs_m", m = m_tisNA )
  if ( !dir.exists( figs_dir ) ) { dir.create( figs_dir, recursive = TRUE ) }
  
  gg_influx <- ggpubr::ggarrange( 
    gg_activ,
    gg_cd69p, 
    
    labels = c( LETTERS[ 1:2 ] ),
    ncol = 2, nrow = 1 )
  pdf_name <- sprintf( "%s/%s_Figure_X_influx.pdf", figs_dir, fig_filename_start )
  ggexport_my( gg_influx, filename = pdf_name, 
               width = 12 * 0.8, height = 6 * 0.8 )

}

for ( celltype in celltypes ) {
  
  # celltype <- "NK"
  # m_tisNA$celltype <- celltype
  # plot_influx( m = m )

}

}
  
  }
