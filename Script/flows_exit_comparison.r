# flows_exit_comparison.r

get_exits_all <- function( m_tisNA = m_tisNA, n_sim = 10000 )
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
    
  
  return( d_exits_all )
}

get_rel_flows_and_exits <- function( m_tisNA = m_tisNA, n_sim = 10000 )
{
  export <- TRUE
  tissues <- get_pmc()$dTissueAllOrderedGroup$tissue.all.ordered[ 1 : 17 ]
  rel_flow_tissues_all <- tibble( celltype = NULL )
  for ( tissue in tissues ) {
    m <- model( m_tisNA, tissue = tissue )
    m$max_lp_chain <- get_max_chain( m )
    chain_dir <- get_path( 
      "chain_dir", m = m, sel_type = "max", sel_chains = m$max_lp_chain )
    al_9 <- get_al_9( m = m )
    Q <- read_Q_matrix( chain_dir )
    if ( ! is.null( Q ) ) {
      rel_flow_tissues <- tibble( celltype = m$celltype )
      rel_flow_tissues$tissue <- tissue
      rel_flow_tissues$al2_Q25 <- al_9[ 2 ] * Q[ 2, 5 ]
      rel_flow_tissues$al2_Q28 <- al_9[ 2 ] * Q[ 2, 8 ]
      rel_flow_tissues$al5_Q52 <- al_9[ 5 ] * Q[ 5, 2 ]
      rel_flow_tissues$al8_Q82 <- al_9[ 8 ] * Q[ 8, 2 ]
      rel_flow_tissues$al3_Q36 <- al_9[ 3 ] * Q[ 3, 6 ]
      rel_flow_tissues$al3_Q39 <- al_9[ 3 ] * Q[ 3, 9 ]
      rel_flow_tissues$al6_Q63 <- al_9[ 6 ] * Q[ 6, 3 ]
      rel_flow_tissues$al9_Q93 <- al_9[ 9 ] * Q[ 9, 3 ]
      
      rel_flow_tissues$al1_Q14 <- al_9[ 1 ] * Q[ 1, 4 ]
      
      rel_flow_tissues$al1_Q12 <- al_9[ 1 ] * Q[ 1, 2 ]
      rel_flow_tissues$al2_Q23 <- al_9[ 2 ] * Q[ 2, 3 ]
      rel_flow_tissues$al3_Q32 <- al_9[ 3 ] * Q[ 3, 2 ]
      
      rel_flow_tissues$al5_Q56 <- al_9[ 5 ] * Q[ 5, 6 ]
      rel_flow_tissues$al6_Q65 <- al_9[ 6 ] * Q[ 6, 5 ]
      
      rel_flow_tissues$al8_Q89 <- al_9[ 8 ] * Q[ 8, 9 ]
      rel_flow_tissues$al9_Q98 <- al_9[ 9 ] * Q[ 9, 8 ]
      
      rel_flow_tissues$al5_Q51 <- al_9[ 5 ] * Q[ 5, 1 ]
      rel_flow_tissues$al8_Q81 <- al_9[ 8 ] * Q[ 8, 1 ]
      rel_flow_tissues$al6_Q61 <- al_9[ 6 ] * Q[ 6, 1 ]
      rel_flow_tissues$al9_Q91 <- al_9[ 9 ] * Q[ 9, 1 ]
      
      rel_flow_tissues$al4_Q41 <- al_9[ 4 ] * Q[ 4, 1 ]
      rel_flow_tissues$al7_Q71 <- al_9[ 7 ] * Q[ 7, 1 ]
      
      rel_flow_tissues_all <- bind_rows( rel_flow_tissues_all, rel_flow_tissues )
    }
  }
  
  d_exits_all <- get_exits_all( m_tisNA = m_tisNA, n_sim = n_sim )
  
  rel_flows_and_exits <- rel_flow_tissues_all %>% 
    mutate( al2_out = al2_Q25 + al2_Q28 + al2_Q23,
            al2_in = al5_Q52 + al8_Q82 + al3_Q32 + al1_Q12,
            al3_out = al3_Q36 + al3_Q39 + al3_Q32,
            al3_in = al6_Q63 + al9_Q93 + al2_Q23 ) %>% 
    left_join( ., d_exits_all, by = c( "tissue" = "tissue" ) ) %>% 
    mutate( die_sum = die_s5 + die_s6 ) %>% 
      dplyr::select( celltype, tissue, al2_out : al3_in, blood_s5 : exit_s4, 
                     die_sum, everything() ) 
    
  
  if ( export ) {
    figs_diag_dir <- get_path( "figs_m_diag", m = m )
    if ( ! dir.exists( figs_diag_dir ) ) { dir.create( figs_diag_dir ) }
    rel_flows_and_exits %>% 
      mutate( across( is.numeric, signif, digits = 4 ) ) %>% 
      write_csv( sprintf( "%s/%s", figs_diag_dir, "rel_flows_and_exits.csv" ) )
  }
  
  return( rel_flows_and_exits )
}
