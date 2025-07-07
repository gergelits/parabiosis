# compare_qij_and_mode.r

if ( FALSE ) {
source( "universal_script_setup.r" )
celltype <- "Tconv"
source( "testing_model_setup.r" )
m$tissue <- "Brain"
}

Q_matrix_to_Q_table <- function( Q )
{
  n <- 6
  Q_table <- expand_grid( i_from = 1 : n, j_to = 1 : n, Q_el = NA ) %>% 
    mutate( i_j = sprintf( "%s_%s", i_from, j_to ) )
  for ( i in 1 : n ) {
    for ( j in 1 : n ) {
      Q_table$Q_el[ Q_table$i_from == i & Q_table$j_to == j ] <- Q[ i, j ]
    }
  }
  return( Q_table )
}



HDI_table_to_Q_table <- function( parabio_fit_HDI )
{
  n <- 6
  Q_table <- expand_grid( i_from = 1 : n, j_to = 1 : n, mode_HDI = NA ) %>% 
    mutate( i_j = sprintf( "%s_%s", i_from, j_to ) )
  for ( i in 1 : n ) {
    for ( j in 1 : n ) {
      s_qij <- sprintf( "q%s%s", i, j )
      s_qij_i1 <- sprintf( "%s_i[1]", s_qij )
      if ( s_qij %in% parabio_fit_HDI$par ) { par_ij <- s_qij } else
        if ( s_qij_i1 %in% parabio_fit_HDI$par ) { par_ij <- s_qij_i1 } else
          par_ij <- NULL
      if ( ! is.null( par_ij ) ) {
        qij_mode <- parabio_fit_HDI %>% 
          filter( par == par_ij ) %>% 
          dplyr::select( "mode" ) %>% 
          unlist() } else
            qij_mode <- 0
    
      Q_table$mode_HDI[ Q_table$i_from == i & Q_table$j_to == j ] <- qij_mode
    }
  }
  return( Q_table )
}


HDI_table_to_Q_matrix <- function( parabio_fit_HDI )
{
  Q_table <- HDI_table_to_Q_table( parabio_fit_HDI )
  n <- max( Q_table$i_from )
  Q_matrix <- matrix( 
    NA, nrow = n, ncol = n )
  for ( i in 1 : n ) {
    for ( j in 1 : n ) {
      Q_matrix[ i, j ] <- 
        Q_table$mode_HDI[ Q_table$i_from == i & Q_table$j_to == j ]
    }
  }
  return( Q_matrix )
}


compare_qij_and_mode <- function( Q, parabio_fit_HDI, qij = "q66" )
{
  Q_table_mat <- Q_matrix_to_Q_table( Q )
  
  if ( ! is.null( qij ) ) {
    i <- substr( qij, 2, 2 )
    j <- substr( qij, 3, 3 )
    cat( "\nQ_matrix.csv\n" )
    print( Q_table_mat %>% filter( i_from == i & j_to == j ) )
    cat( "\nparabio_fit_HDI\n" )
    print( parabio_fit_HDI %>% 
             filter( par == sprintf( "%s_i[1]", qij ) ) %>% 
             dplyr::select( "mode", "mean" ) )
  }
  
  if ( is.null( qij ) ) {
    Q_table_HDI <- HDI_table_to_Q_table( parabio_fit_HDI )
    Q_tables_diff <- Q_table_mat %>% 
      dplyr::select( -i_from, -j_to ) %>% 
      left_join( ., Q_table_HDI %>% dplyr::select( -i_from, -j_to ), 
                 by = c( "i_j" = "i_j" ) ) %>% 
      mutate( qij_diff = Q_el - mode_HDI,
              qij_rat = Q_el / mode_HDI )                                         
    return( Q_tables_diff )
  }
}


if ( TESTING <- FALSE ) {
  source( "universal_script_setup.r" )
  celltype <- "Tconv"
  source( "testing_model_setup.r" )
  m$tissue <- "Brain"  
model_dir_r <- get_path( "model_dir_r", m = m )
chain_dir <- get_path( "chain_dir", m = m, sel_type = "max", sel_chains = "1" )
Q <- read_Q_matrix( chain_dir )
parabio_fit_HDI <- read_csv( 
  sprintf( "%s/parabio_fit_HDI.csv",
           get_path( "model_dir", m = m, model_dir_r = model_dir_r ) ) )

compare_qij_and_mode( Q = Q, parabio_fit_HDI = parabio_fit_HDI, qij = NULL ) %>% 
  arrange( qij_rat ) %>% 
  as.data.frame()

compare_qij_and_mode( Q = Q, parabio_fit_HDI = parabio_fit_HDI, qij = "q61" )
compare_qij_and_mode( Q = Q, parabio_fit_HDI = parabio_fit_HDI, qij = "q63" )
compare_qij_and_mode( Q = Q, parabio_fit_HDI = parabio_fit_HDI, qij = "q65" )
compare_qij_and_mode( Q = Q, parabio_fit_HDI = parabio_fit_HDI, qij = "q66" )

}
