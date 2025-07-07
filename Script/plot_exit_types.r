# plot_exit_types.r

plot_exit_types <- function( m_tisNA, n_sim = 10000 )
{
  tissues <- get_pmc()$dTissueAllOrderedGroup$tissue.all.ordered[ 1 : 17 ]
  d_exits_all_tmp <- tibble( tissue = NULL )
  for ( tissue in tissues ) {
    m <- model( m_tisNA, tissue = tissue )
    fig_fname_start <- get_fig_fname_start( m = m )
    figs_dwell_dir <- get_path( "figs_m_dwell", m = m )
    prop_exit_file <- sprintf( "%s/%s_%s_%isims_prop_exit.csv", figs_dwell_dir, 
                               fig_fname_start, tolower( m$tissue ), n_sim )
    if ( file.exists( prop_exit_file ) ) {
      d_exits <- read_csv( prop_exit_file ) %>% 
        mutate( tissue = m$tissue ) %>% 
        dplyr::select( tissue, state, freq ) %>% 
        pivot_wider( names_from = state, values_from = freq ) 
      d_exits_all_tmp <- bind_rows( d_exits_all_tmp, d_exits )
    }
  }
  
  d_exits_all <- 
    d_exits_all_tmp %>% 
    mutate( die_s5 = ifelse( is.na( die_s5 ), 0, die_s5 ), 
            die_s6 = ifelse( is.na( die_s6 ), 0, die_s6 ), 
            exit_s4 = ifelse( is.na( exit_s4 ), 0, exit_s4 ) )
    
  d_exits_groups <- 
    d_exits_all %>% 
    pivot_longer( cols = c( "blood_s5", "blood_s6", "die_s5", "die_s6", 
                            "exit_s4" ),
                  names_to = "exit_type", values_to = "proportion" ) %>% 
    mutate( exit_type = ifelse( exit_type == "blood_s5", 
                                "to_blood_activ", exit_type ),
            exit_type = ifelse( exit_type == "blood_s6", 
                                "to_blood_cd69p", exit_type ),
            exit_type = ifelse( exit_type == "die_s5", "die_activ", exit_type ),
            exit_type = ifelse( exit_type == "die_s6", "die_cd69p", exit_type ),
            exit_type = ifelse( exit_type == "exit_s4", 
                                "to_blood_or_die_naive", exit_type ) ) %>% 
    
    mutate( exit_type = as.factor( exit_type ) ) %>% 
    left_join( ., get_pmc()$dTissueAllOrderedGroup.f %>% 
                 dplyr::select( tissue.all, f.tissue, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) )
    
    
  xintercept1 <- d_exits_groups$f.tissue[ 
    d_exits_groups$f.tissue.group == "BoneMarrow" ] %>% unique() %>% length()
  xintercept2 <- d_exits_groups$f.tissue[ 
    d_exits_groups$f.tissue.group == "Lymphoid" ] %>% unique() %>% length()
  xintercept3 <- d_exits_groups$f.tissue[ 
    d_exits_groups$f.tissue.group == "Non-lymphoid" ] %>% unique() %>% length()
  
  cm_to_in <- 1 / 2.54   
  f <- cm_to_in * 1.7
  
  if ( m$celltype != "Treg" ) {
    exit_labels <- c( 
      "Activated (die)", "CD69+ (die)", 
      "Activated (exit)", "CD69+ (exit)", "Naive (die/exit)" )
    } else {
      exit_labels <- c( 
        "Antigen-experienced (die)", "CD69+ (die)", "Antigen-experienced (exit)", 
        "CD69+ (exit)", "Resting (die/exit)" ) 
      }
  
  gg_exit_types <- 
    d_exits_groups %>% 
    ggplot( aes( x = f.tissue, y = proportion, fill = exit_type ) ) +
    geom_bar( stat = "identity", position = "dodge" ) +
    geom_vline( xintercept = c( xintercept1 + 0.5, 
                                xintercept1 + xintercept2 + 0.5,
                                xintercept1 + xintercept2 + xintercept3 + 0.5 ),
                linetype = "dashed" ) +
    labs( y = "Proportion", fill = "Exit type" ) +
    scale_fill_discrete( labels = exit_labels ) +
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
           strip.text.x = element_text( size = 12 * f ) ); gg_exit_types
    
    
    fig_fname_start <- get_fig_fname_start( m = m_tisNA ) 
    figs_dwell_dir <- get_path( "figs_m", m = m_tisNA )
    pdf_name <- sprintf( "%s/%s_%isims_exit_types.pdf", 
                         figs_dwell_dir, fig_fname_start, n_sim )
    
    ggexport_my( gg_exit_types, filename = pdf_name, 
                 width = 29.7 * cm_to_in, height = 10 * cm_to_in )
}
