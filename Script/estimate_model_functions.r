# estimate_model_functions.r

read_Q_matrix <- function( chain_dir )
{
  Q_matrix_files <- list.files( 
    chain_dir, pattern = "^Q_matrix.*.csv", full.names = TRUE )
  if ( length( Q_matrix_files ) == 1 ) {
    Q <- read.csv( Q_matrix_files[ 1 ] )[ , -1 ] %>% as.matrix
    rownames( Q ) <- colnames( Q ) } else {
      Q <- NULL
    }
  return( Q )
}


get_parabio_fit_HDI <- function( m, chain = NA )
{
  chain <- ifelse( ! is.na( chain ), chain, m$chain )
  chain_dir <- get_path( "chain_dir", m = m, chain = chain,
                         sel_type = "each", sel_chains = chain )
  
  parabio_fit_HDI_f <- file.path( chain_dir, "parabio_fit_HDI.csv" )
  if ( ! file.exists( parabio_fit_HDI_f ) )
  {
    calc_HDI <- calculate_HDI( m = m, chain = chain )
    if ( is.null( calc_HDI ) ) { return( NULL ) }
  }
  parabio_fit_HDI <- read_csv( parabio_fit_HDI_f )
  return( parabio_fit_HDI )
}


read_af_matrix <- function( m, chain = NA )
{
  chain <- ifelse( ! is.na( chain ), chain, m$chain )
  parabio_fit_HDI <- get_parabio_fit_HDI( m = m, chain = chain )
  af_tmp <- parabio_fit_HDI %>% 
    dplyr::select( par, mode ) %>% 
    filter( row_number() %in% grep( pattern = "^af\\[", x = par ) ) %>% 
    mutate( af_ind_1 = gsub( "af\\[(.*)\\,(.*)\\]", "\\1", par ) ) %>% 
    mutate( af_ind_2 = gsub( "af\\[(.*)\\,(.*)\\]", "\\2", par ) ) %>% 
    arrange( af_ind_2, af_ind_1 )
  af_matrix <- matrix( af_tmp$mode, 
                       nrow = af_tmp$af_ind_1 %>% as.integer %>% max, 
                       ncol = af_tmp$af_ind_2 %>% as.integer %>% max )
  return( af_matrix )
}


plot_par_densities_and_calc_Q <- function( 
    m, parabio.fit, sel_type, sel_chains, iter = "all" )
{       
  ps <- m$ps
  model_name <- m$model_name; 
  model_dir <- get_path( "model_dir", m = m )
  model_vars_dir <- get_path( "model_vars_dir", m = m )
  
  chain_dir <- get_path( 
    "chain_dir", m = m, model_dir = model_dir, 
    sel_type = sel_type, sel_chains = sel_chains )
    
  if ( ! dir.exists( chain_dir ) ) { 
    dir.create( chain_dir, recursive = TRUE ) }
  
  # pb.tissue, pb.hostdonor.popul, parabio.data, pb.hostdonor.tissue.popul, 
  # pb.popul
  load( file = sprintf( "%s/pb.RData", model_vars_dir ) )
  # other_vars = parabio.data
  load( file = sprintf( "%s/other_vars.RData", model_vars_dir ) )
  # other_vars needed for tn_
  load( file = sprintf( "%s/model_vars_.RData", model_vars_dir ) )
  
  if ( GET_FITTED_PARAMETERS_FROM_MCMC_SAMPLE_TO_R_ENV <- TRUE ) {
  
    if( ! is.character( sel_chains ) ) { browser() }
    
      sel_chains <- stringr::str_split_1( sel_chains, pattern = "_" ) %>% 
        as.numeric() 
    # }
    # else: sel_chains is a single integer
    
    if ( iter == "all" ) {
      # get MCMC sample    
      parabio.fit_draws <- posterior::as_draws( parabio.fit )
      draws_ok <- posterior::subset_draws( 
        parabio.fit_draws, chain = sel_chains )
      summary_draws_ok <- posterior::summarise_draws( draws_ok )
      write_csv( summary_draws_ok, 
                 file = sprintf( "%s/parabio_fit_print_nOkCh=%s.csv", 
                                 chain_dir, sum( !is.na( sel_chains ) ) ) )
    }
    parabio.fit.sample <- rstan::extract( parabio.fit, permuted = FALSE )
    
    if ( GET_ALL_THE_CALCULATED_PARAMS_TO_THE_ENV <- TRUE ) {
      # get the fitted parameters ("free" parameters here only) from the MCMC 
      # sample for sel_chains
      
      if ( ASSIGN_VALUES_TO_FREE_PARAMS__BACKTICKS_R <- TRUE ) {
        for ( mp in dimnames( parabio.fit.sample )$parameters )
        {
          if ( iter == "all" ) {
            pf.sample <- as.vector( parabio.fit.sample[ , sel_chains, mp ] )
            if ( any( is.na( pf.sample ) ) ) { browser( ) }
            
            pf.sample.mode <- mlv( pf.sample, method = "venter", type = "shorth" )
            pf.sample.density <- density( pf.sample, n = 1000, cut = 0 )
            pf.sample.density.max.x <- pf.sample.density$x[
              which.max( pf.sample.density$y ) ]
            assign( mp, pf.sample.density.max.x )
            
            
            png( filename = sprintf( "%s/density_%s.png", chain_dir, mp ),
                 width = 1024, height = 1024 )
            par( mfrow = c( 2, 1 ), mar = c( 5, 6, 1, 1 ), lwd = 3, cex.lab = 2,
                 cex.axis = 2 )
            hist( pf.sample, breaks = 30, freq = FALSE, xlab = mp, main = "" )
            abline( v = pf.sample.density.max.x, col = "blue3" )
            abline( v = pf.sample.mode, col = "red3", lty = 2 )
            plot( pf.sample.density, xlab = mp, main = "" )
            abline( v = pf.sample.density.max.x, col = "blue3" )
            abline( v = pf.sample.mode, col = "red3", lty = 2 )
            dev.off()
          } else {
            stopifnot( iter <= dim( parabio.fit.sample )[ 1 ] )
            pf.value.iter <- as.vector( parabio.fit.sample[ , sel_chains, mp ] )[ 
              as.numeric( iter ) ]
            assign( mp, pf.value.iter )
          }
        }
      }
      
      # parse the .stan code to create the whole Q matrix in the R environment
      stan.file <- readr::read_lines( sprintf( "%s/%s.stan", ps$CODE_PATH, model_name ) )
      
      if ( FIND_STAN_CODE_AND_INITIALIZE_VECTORS__STANDARD_R <- TRUE ) {
        # get the .stan code part, which will initialize vectors in R
        # ( in .stan it is done by different code in part parameters{} )
        
        # get the .stan code part
        n.line.definitions.start <- stan.file %>% str_which(
          "\\/\\/ start of definitions for R", negate = FALSE )
        n.line.definitions.end <- stan.file %>% str_which(
          "\\/\\/ end of definitions for R", negate = FALSE )
        
        # initialize vectors
        for ( n.line in ( n.line.definitions.start + 1 ) : ( n.line.definitions.end - 1 ) ) {
          n.line.expression <- 
            stan.file[ n.line ] %>% 
            str_replace_all( ., pattern = "\\{\\{//REMOVE IN R", 
                             replacement = c( "{{//REMOVE IN R" = "" ) ) %>%         # replace .stan min ( {a, b} ) syntax
            str_replace_all( ., pattern = "\\}\\}//REMOVE IN R", 
                             replacement = c( "}}//REMOVE IN R" = "" ) ) %>%         # replace .stan min ( {a, b} ) syntax
            
            gsub( "\\/\\/\\#" , "\\#", . ) %>%                                       # remove stan comments but keep R comments
            gsub( "\\/\\/" , "", . ) %>%                                             # replace stan comments - hack to "send" some extra R code to the R parser
            
            str_replace_all( ., pattern = "\\{\\{", 
                             replacement = c( "{{" = "XXTMP1" ) ) %>%                # replace .stan min ( {a, b} ) syntax
            str_replace_all( ., pattern = "\\}\\}", 
                             replacement = c( "}}" = "XXTMP2" ) ) %>%                # replace .stan min ( {a, b} ) syntax
            
            str_replace_all( ., pattern = "\\{", replacement = c( "{" = "c(" ) ) %>% # replace .stan min ( {a, b} ) syntax
            str_replace_all( ., pattern = "\\}", replacement = c( "}" = ")" ) ) %>%  # replace .stan min ( {a, b} ) syntax
            str_replace_all( ., pattern = "XXTMP1", 
                             replacement = c( "XXTMP1" = "{" ) ) %>%                 # replace .stan min ( {a, b} ) syntax
            str_replace_all( ., pattern = "XXTMP2", 
                             replacement = c( "XXTMP1" = "}" ) ) %>%                 # replace .stan min ( {a, b} ) syntax
            
            str_replace_all( ., pattern = "int", replacement = c( "int" = "" ) ) %>%
            str_replace_all( ., pattern = "real", replacement = c( "real" = "" ) ) %>%
            
            parse( text = . )
            
          eval( n.line.expression )
        } 
      }
      
      # Assign the modes of !!VECTORED & SAMPLED!! (ONLY???) values from .stan model to R variables
      
      if ( CONVERT_BACKTICKS_R_TO_STANDARD_R <- TRUE ) {
        
        if ( MODEL_1000p <- TRUE ) {
          for ( i in 1 : ( tn_-1 ) ) {
            if ( exists( "q14_simplex_i[1]" ) ) {                                      
              q14_simplex_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q14_simplex_i[", i, "]`" ) ) %>% eval() }
          }
          
          for ( i in 1 : ( tn_-1 ) ) {
            q51_i[ i ] <- parse( text = sprintf( "%s%i%s", "`q51_i[", i, "]`" ) ) %>% eval()
            q52_i[ i ] <- parse( text = sprintf( "%s%i%s", "`q52_i[", i, "]`" ) ) %>% eval()
            q41_i[ i ] <- parse( text = sprintf( "%s%i%s", "`q41_i[", i, "]`" ) ) %>% eval()
            
            # since m0050
            if ( exists( "q25_min_N5_i[1]" ) ) {                                      
              q25_min_N5_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q25_min_N5_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q25_max_N5_i[1]" ) ) {                                      
              q25_max_N5_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q25_max_N5_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q52_min_N5_i[1]" ) ) {                                      
              q52_min_N5_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q52_min_N5_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q52_max_N5_i[1]" ) ) {                                      
              q52_max_N5_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q52_max_N5_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q52_min_N2_i[1]" ) ) {                                      
              q52_min_N2_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q52_min_N2_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q52_max_N2_i[1]" ) ) {                                      
              q52_max_N2_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q52_max_N2_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q52_min_i[1]" ) ) {                                      
              q52_min_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q52_min_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q52_max_i[1]" ) ) {                                      
              q52_max_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q52_max_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q52_base_i[1]" ) ) {                                      
              q52_base_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q52_base_i[", i, "]`" ) ) %>% eval() }

                        
            if ( exists( "q61_q65_sum_i[1]" ) ) {                                      
              q61_q65_sum_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q61_q65_sum_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q61_q65_sum_propos_i[1]" ) ) {                                      
              q61_q65_sum_propos_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q61_q65_sum_propos_i[", i, "]`" ) ) %>% eval() }
            
            
            if ( exists( "q14_min_I_i[1]" ) ) {                                      
              q14_min_I_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q14_min_I_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q14_max_I_i[1]" ) ) {                                      
              q14_max_I_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q14_max_I_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q14_base_i[1]" ) ) {                                      
              q14_base_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q14_base_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q14_propos_i[1]" ) ) {                                      
              q14_propos_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q14_propos_i[", i, "]`" ) ) %>% eval() }
            
            
            
            
            if ( exists( "q14_max_i[1]" ) ) {                                      
              q14_max_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q14_max_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q41_u_code_i[1]" ) ) {                                      
              q41_u_code_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q41_u_code_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q45_u_code_i[1]" ) ) {                                      
              q45_u_code_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q45_u_code_i[", i, "]`" ) ) %>% eval() }
            
            
            if ( exists( "q36_propos_i[1]" ) ) {                                      
              q36_propos_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q36_propos_i[", i, "]`" ) ) %>% eval() }
            
            if ( exists( "q51_u_code_i[1]" ) ) {                                      
              q51_u_code_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q51_u_code_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q61_u_code_i[1]" ) ) {                                      
              q61_u_code_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q61_u_code_i[", i, "]`" ) ) %>% eval() }
            
            if ( exists( "al1_q14_min_i[1]" ) ) {                                      
              al1_q14_min_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`al1_q14_min_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "al4_q41_flow_i[1]" ) ) {                                      
              al4_q41_flow_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`al4_q41_flow_i[", i, "]`" ) ) %>% eval() }
            
            
            if ( exists( "q41_flow_i[1]" ) ) {                                      
              q41_flow_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q41_flow_i[", i, "]`" ) ) %>% eval() }
            
            if ( exists( "wgh4_free_i[1]" ) ) {                                      
              wgh4_free_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`wgh4_free_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "wgh5_free_i[1]" ) ) {                                      
              wgh5_free_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`wgh5_free_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "wgh6_free_i[1]" ) ) {                                      
              wgh6_free_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`wgh6_free_i[", i, "]`" ) ) %>% eval() }
            
            if ( exists( "wgh4_i[1]" ) ) {                                      
              wgh4_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`wgh4_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "wgh5_i[1]" ) ) {                                      
              wgh5_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`wgh5_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "wgh6_i[1]" ) ) {                                      
              wgh6_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`wgh6_i[", i, "]`" ) ) %>% eval() }
            
            if ( exists( "q14_min_i[1]" ) ) {                                      
              q14_min_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q14_min_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q14_i[1]" ) ) {                                      
              q14_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q14_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q25_propos_i[1]" ) ) {                                      
              q25_propos_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q25_propos_i[", i, "]`" ) ) %>% eval() }
            
            
            if ( exists( "q41_base_i[1]" ) ) {                                      
              q41_base_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q41_base_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q41a_i[1]" ) ) {                                      
              q41a_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q41a_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q41b_i[1]" ) ) {                                      
              q41b_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q41b_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q36_sum_max_simplex_i[1]" ) ) {                                      
              q36_sum_max_simplex_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q36_sum_max_simplex_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q36_sum_min_simplex_i[1]" ) ) {                                      
              q36_sum_min_simplex_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q36_sum_min_simplex_i[", i, "]`" ) ) %>% eval() }
            
            if ( exists( "q36_max_S_i[1]" ) ) {                                      
              q36_max_S_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q36_max_S_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q36_max_I_i[1]" ) ) {                                      
              q36_max_I_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q36_max_I_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q36_max_i[1]" ) ) {                                      
              q36_max_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q36_max_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q36_min_S_i[1]" ) ) {                                      
              q36_min_S_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q36_min_S_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q36_min_I_i[1]" ) ) {                                      
              q36_min_I_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q36_min_I_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q36_min_i[1]" ) ) {                                      
              q36_min_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q36_min_i[", i, "]`" ) ) %>% eval() }
            
            if ( exists( "q36_base_i[1]" ) ) {                                              
              q36_base_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q36_base_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q36_simplex_i[1]" ) ) {                                      
              q36_simplex_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q36_simplex_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q36_i[1]" ) ) {                                      
              q36_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q36_i[", i, "]`" ) ) %>% eval() }
            
            
            
            if ( exists( "q61_i[1]" ) ) {                                              
              q61_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q61_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q63_i[1]" ) ) {                                              
              q63_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q63_i[", i, "]`" ) ) %>% eval() }
            
            if ( exists( "q65_base_i[1]" ) ) {                                         
              q65_base_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q65_base_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q65_i[1]" ) ) {                                              
              q65_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q65_i[", i, "]`" ) ) %>% eval() }
            
            if ( exists( "q56_base_i[1]" ) ) {                                         
              q56_base_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q56_base_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q56_i[1]" ) ) {                                         
              q56_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q56_i[", i, "]`" ) ) %>% eval() }
            
            
            if ( exists( "q25_S_i[1]" ) ) {                                      
              q25_S_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q25_S_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q25_inminmax_simplex_i[1]" ) ) {                                      
              q25_inminmax_simplex_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q25_inminmax_simplex_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "al_2_q25__al_1[1]" ) ) {                                      
              al_2_q25__al_1[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`al_2_q25__al_1[", i, "]`" ) ) %>% eval() }
            
            
            if ( exists( "q51_min_i[1]" ) ) {                                      
              q51_min_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q51_min_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q51_max_i[1]" ) ) {                                      
              q51_max_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q51_max_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q51_base_i[1]" ) ) {                                      
              q51_base_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q51_base_i[", i, "]`" ) ) %>% eval() }
            
            
            if ( exists( "q45_max_i[1]" ) ) {                                      
              q45_max_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q45_max_i[", i, "]`" ) ) %>% eval() }
            
            
            if ( exists( "q25_sum_max_simplex_i[1]" ) ) {                                      
              q25_sum_max_simplex_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q25_sum_max_simplex_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q25_sum_min_simplex_i[1]" ) ) {                                      
              q25_sum_min_simplex_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q25_sum_min_simplex_i[", i, "]`" ) ) %>% eval() }
            
            if ( exists( "q25_max_S_i[1]" ) ) {                                      
              q25_max_S_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q25_max_S_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q25_max_I_i[1]" ) ) {                                      
              q25_max_I_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q25_max_I_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q25_max_i[1]" ) ) {                                      
              q25_max_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q25_max_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q25_min_S_i[1]" ) ) {                                      
              q25_min_S_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q25_min_S_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q25_min_I_i[1]" ) ) {                                      
              q25_min_I_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q25_min_I_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q25_min_i[1]" ) ) {                                      
              q25_min_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q25_min_i[", i, "]`" ) ) %>% eval() }
            
            if ( exists( "q25_base_i[1]" ) ) {                                              
              q25_base_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q25_base_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q25_simplex_i[1]" ) ) {                                      
              q25_simplex_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q25_simplex_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q25_i[1]" ) ) {                                      
              q25_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q25_i[", i, "]`" ) ) %>% eval() }
            
            if ( exists( "q54_base_i[1]" ) ) {                                      
              q54_base_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q54_base_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q54_max_i[1]" ) ) {                                      
              q54_max_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q54_max_i[", i, "]`" ) ) %>% eval() }
            if ( exists( "q54_i[1]" ) ) {                                      
              q54_i[ i ] <- 
                parse( text = sprintf( "%s%i%s", "`q54_i[", i, "]`" ) ) %>% eval() }
            
          }  
          if ( exists( "q12" ) ) {                                                     
            q12 <- parse( text = sprintf( "`q12`" ) ) %>% eval() 
          } 
        }
      }
      
      if ( FIND_STAN_CODE_AND_ASSIGN_VALUES_TO_NONFREE_PARAMETERS__STANDARD_R <- TRUE ) {
        # get the .stan code part, which will assign values to the non-free parameters
        n.line.nonfree.start <- stan.file %>% str_which( 
          "\\/\\/ start of nonfree variables", negate = FALSE )
        n.line.nonfree.end <- stan.file %>% str_which( 
          "\\/\\/ end of nonfree variables", negate = FALSE )
        
        # assign the value to the non-free parameters based on the extracted .stan code
        for ( n.line in ( n.line.nonfree.start + 1 ) : ( n.line.nonfree.end - 1 ) ) {
          n.line.expression <- 
            stan.file[ n.line ] %>%
            str_replace_all( ., pattern = "\\{\\{//REMOVE IN R", 
                             replacement = c( "{{//REMOVE IN R" = "" ) ) %>%         # replace .stan min ( {a, b} ) syntax
            str_replace_all( ., pattern = "\\}\\}//REMOVE IN R", 
                             replacement = c( "}}//REMOVE IN R" = "" ) ) %>%         # replace .stan min ( {a, b} ) syntax
            
            gsub( "\\/\\/\\/\\/" , "#", . ) %>%                                    # CORRECTION NEEDED due to not allowing # in .stan
            gsub( "\\/\\/\\#" , "\\#", . ) %>%                                       # remove stan comments but keep R comments
            gsub( "\\/\\/" , "", . ) %>%                                             # replace stan comments - hack to "send" some extra R code to the R parser
            
            str_replace_all( ., pattern = "\\{\\{", 
                             replacement = c( "{{" = "XXTMP1" ) ) %>%                # replace .stan min ( {a, b} ) syntax
            str_replace_all( ., pattern = "\\}\\}", 
                             replacement = c( "}}" = "XXTMP2" ) ) %>%                # replace .stan min ( {a, b} ) syntax
            
            
            str_replace_all( ., pattern = "\\{", replacement = c( "{" = "c(" ) ) %>% # replace .stan min ( {a, b} ) syntax
            str_replace_all( ., pattern = "\\}", replacement = c( "}" = ")" ) ) %>%  # replace .stan min ( {a, b} ) syntax
            str_replace_all( ., pattern = "XXTMP1", 
                             replacement = c( "XXTMP1" = "{" ) ) %>%                 # replace .stan min ( {a, b} ) syntax
            str_replace_all( ., pattern = "XXTMP2", 
                             replacement = c( "XXTMP1" = "}" ) ) %>%                 # replace .stan min ( {a, b} ) syntax
            
            str_replace_all( ., pattern = "int", replacement = c( "int" = "" ) ) %>%
            str_replace_all( ., pattern = "real", replacement = c( "real" = "" ) ) %>%
            
            parse( text = . )
            
          eval( n.line.expression )
        } 
      }
    }
    
  }
  
  # having variables in environment build intensity matrix Q
  if ( HAVING_PARAMS_IN_STANDARD_R_ENV_BUILD_AND_SAVE_INTENSITY_MATRIX_Q <- TRUE ) {
    Q <- matrix( 
      0, nrow = 2 * pb.tissue.n * pb.popul.n, ncol = 2 * pb.tissue.n * pb.popul.n )
    
    n.line.Q.alignment.start <- stan.file %>% 
      str_which( ., "\\/\\/ start of Q alignment", negate = FALSE )
    n.line.Q.alignment.end <- stan.file %>% 
      str_which( ., "\\/\\/ end of Q alignment", negate = FALSE )
    stan.file[ ( n.line.Q.alignment.start + 1 ) : ( n.line.Q.alignment.end - 1 ) ]
    
    for ( n.line in ( n.line.Q.alignment.start + 1 ) : ( n.line.Q.alignment.end - 1 ) ) {
      n.line.expression <- stan.file[ n.line ] %>% 
        # replace stan comments
        gsub( "\\/\\/" , "\\#", . ) %>%
        # replace .stan min ( {a, b} ) syntax
        str_replace_all( ., pattern = "\\{\\{", 
                         replacement = c( "{{" = "XXTMP1" ) ) %>%                
        # replace .stan min ( {a, b} ) syntax
        str_replace_all( ., pattern = "\\}\\}", 
                         replacement = c( "}}" = "XXTMP2" ) ) %>%                
        # replace stan min ( {a, b} ) syntax
        str_replace_all( ., pattern = "\\{", replacement = c( "{" = "c(" ) ) %>% 
        # replace stan min ( {a, b} ) syntax
        str_replace_all( ., pattern = "\\}", replacement = c( "}" = ")" ) ) %>%  
        # replace stan vector [a, b] syntax
        str_replace_all( ., pattern = "= \\[", 
                         replacement = c( "= \\[" = "= c(" ) ) %>%               
        str_replace_all( ., pattern = "\\];", 
                         replacement = c( "\\];" = ");" ) ) %>%
        # replace .stan min ( {a, b} ) syntax
        str_replace_all( ., pattern = "XXTMP1", 
                         replacement = c( "XXTMP1" = "{" ) ) %>%                 
        # replace .stan min ( {a, b} ) syntax
        str_replace_all( ., pattern = "XXTMP2", 
                         replacement = c( "XXTMP1" = "}" ) ) %>%                 
        parse( text = . )
        
      eval( n.line.expression )
    } 
    
    
    dimnames( Q ) <- list( pb.hostdonor.tissue.popul, pb.hostdonor.tissue.popul )
    if ( iter == "all" ) {
      error.suffix <- ""
      if ( ! ( all ( Q - diag( diag( Q ) ) >= 0 ) ) ) {
        error.suffix <- paste0( error.suffix, "_Qsigns" ) }
      if ( ! ( all ( abs( rowSums( Q ) ) < 1e-10 ) ) ) {
        error.suffix <- paste0( error.suffix, "_rowSums0" ) }
      write.csv( Q, file = sprintf( 
        "%s/Q_matrix%s.csv", chain_dir, error.suffix ), quote = FALSE )
    } else {
      write.csv( Q, file = sprintf( 
        "%s/tmp_Q_matrix%s.csv", chain_dir, "_iter_x" ), quote = FALSE )
    }
  }
  return( Q )
}  


get_and_plot_trajectories <- function( 
  m = NULL,
  sel_chains, sel_type, 
  sel_models_file_name )
{
  celltype <- m$celltype; ps <- m$ps; mcmc_pars <- m$mcmc_pars
  model_dir <- get_path( "model_dir", m = m )
  model_vars_dir <- get_path( "model_vars_dir", m = m )
  
  load( file = sprintf( "%s/pb.RData", model_vars_dir ) )
  load( file = sprintf( "%s/other_vars.RData", model_vars_dir ) )              # other_vars = parabio.data
  # loading for w0i_ etc.
  load( file = sprintf( "%s/model_vars_.RData", model_vars_dir ) )
  chain_dir <- get_path( "chain_dir", m = m, model_dir = model_dir, 
                         sel_type = sel_type, sel_chains = sel_chains )
  Q <- read_Q_matrix( chain_dir )
  
  
  # get limiting distribution in paired mice
  if ( FROM_Q_GET_LIMITING_DISTRIBUTION_FOR_TRAJECTORIES_PLOT <- TRUE ) {
    Qs <- Q
    n.zeroed <- 0
    Qs[ , 1 ] <- rep( 1, 2 * pb.tissue.n * pb.popul.n - n.zeroed )
    bl <- solve( t( Qs ), c( 1, rep( 0, 2 * pb.tissue.n * pb.popul.n - 1 - n.zeroed ) ) )
    names( bl )[ c( 3, pb.tissue.n * pb.popul.n + 3 ) ] <- 
      rownames( Q )[ c( 3, pb.tissue.n * pb.popul.n + 3 ) ]
  }
  
  # get average trajectories
  # INPUT (specific for sel_chains): Q
  if ( GET_AVERAGE_TRAJECTORIES <- TRUE ) {
    week.calc <- seq( w0i_, pb.figure.week.max, length.out = 200 )
    alpha.calc <- t( sapply( week.calc, function( w )
      as.vector( a0i_ %*% expm( 7 * ( w - w0i_ ) * Q ) ) ) )
    y.calc <- lapply( pb.tissue, function( tis ) {
      tissue.idx <- grep( sprintf( "\\.%s\\.", tis ), pb.hostdonor.tissue.popul )
      yc <- alpha.calc[ , tissue.idx ]
      yc <- sweep( yc, 1, rowSums( yc ), "/" )
      colnames( yc ) <- pb.hostdonor.popul
      yc
    } )
    names( y.calc ) <- pb.tissue
  }
  
  # plot trajectories
  if ( PLOT_TRAJECTORIES <- TRUE ) {
    l_ggfigs <- list()
    l_ggfigs_ext <- list()
    
    for ( tis in pb.tissue )
    { 
      l_ggfigs[[ which( pb.tissue == tis ) ]] <- list()
      l_ggfigs_ext[[ which( pb.tissue == tis ) ]] <- list()
      names( parabio.data ) <- gsub( 
        pattern = "celltype", replacement = tolower( celltype ), 
        names( parabio.data ) )
      tissue.week <- parabio.data[ parabio.data$tissue == tis, "week" ]
      tissue.data <- parabio.data[ parabio.data$tissue == tis,
                                   pb.hostdonor.popul ]
      tissue.data <- sweep( tissue.data, 1, rowSums( tissue.data ), "/" )
      
      limit.data <- bl[ grep( sprintf( "\\.%s\\.", tis ),
                              pb.hostdonor.tissue.popul ) ]
      limit.data <- limit.data / sum( limit.data )
      names( limit.data ) <- pb.hostdonor.popul
      
      for ( pop in pb.popul )
      {
        host.pop <- sprintf( "host.%s", pop )
        donor.pop <- sprintf( "donor.%s", pop )
        ggdata.calc <- data.frame( x = week.calc,
                                   y1 = y.calc[[ tis ]][ , host.pop ],
                                   y2 = y.calc[[ tis ]][ , donor.pop ] )
        ggdata.exp <- data.frame( x = tissue.week,
                                  y1 = tissue.data[[ host.pop ]], 
                                  y2 = tissue.data[[ donor.pop ]] )
        ggdata.exp.mean <- data.frame( t( sapply( pb.week, function( w )
          c( w, mean( tissue.data[ tissue.week == w, host.pop ] ),
             mean( tissue.data[ tissue.week == w, donor.pop ] ) ) ) ) )
        names( ggdata.exp.mean ) <- c( "x", "y1", "y2" )
        # y.max <- 1
        
        ggfig_pre <- ggplot( ggdata.calc ) +
          geom_hline( yintercept = limit.data[ host.pop ], linetype = "dashed" ) +
          geom_hline( yintercept = limit.data[ donor.pop ], linetype = "dotted" ) +
          geom_line( aes( x, y1 ), color = "blue3" ) +
          geom_line( aes( x, y2 ), color = "red3" ) +
          geom_point( aes( x, y1 ), data = ggdata.exp, size = 1,
                      position = position_nudge( - 0.1 ), color = "blue3" ) +
          geom_point( aes( x, y2 ), data = ggdata.exp, size = 1,
                      position = position_nudge( 0.1 ), color = "red3" ) +
          geom_point( aes( x, y1 ), data = ggdata.exp.mean, size = 3,
                      position = position_nudge( - 0.1 ), color = "blue3" ) +
          geom_point( aes( x, y2 ), data = ggdata.exp.mean, size = 3,
                      position = position_nudge( 0.1 ), color = "red3" ) +
          scale_x_continuous( breaks = seq( 0, pb.figure.week.max, by = 2 ),
                              name = "Week" ) +
          ggtitle( sprintf( 
            "%s %s", 
            ifelse( pb.tissue.label[ tis ] == "BoneMarrow", 
                    "Bone marrow", pb.tissue.label[ tis ] ),
            ifelse( ( pb.popul.label[ pop ] == "Treg Naive" ), 
                    "Treg Resting", pb.popul.label[ pop ] ) ) ) +           
          theme_bw() +
          theme( panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.text = element_text( size = 14 ),
                 axis.title = element_text( size = 14 ),
                 plot.title = element_text( size = 16 ) )
        
        for ( y.max in c( 1, NA ) )
        {
          ggfig <- ggfig_pre +
            scale_y_continuous( 
              limits = c( 0, y.max ), 
              name = sprintf( 
                "Fraction in %s %s", 
                ifelse( pb.tissue.label[ tis ] == "BoneMarrow", 
                        "Bone marrow", pb.tissue.label[ tis ] ), 
                pb.parent.popul.label ) )
          if ( is.na( y.max ) ) {
            ggsave( sprintf( "%s/%s_%s_ext.png", chain_dir, tis, pop ), ggfig, 
                    width = 8, height = 6 )
            l_ggfigs_ext[[ which( pb.tissue == tis ) ]][[ 
              which( pb.popul == pop ) ]] <- ggfig
          } else if ( y.max == 1 ) {
            ggsave( sprintf( "%s/%s_%s.png", chain_dir, tis, pop ), ggfig, 
                    width = 8, height = 6 )
            l_ggfigs[[ which( pb.tissue == tis ) ]][[ 
              which( pb.popul == pop ) ]] <- ggfig
          }
        }
      }
    }
    
    # all_trajectories figure for model of tissue_i
    if ( length( pb.tissue ) == 4 ) {
      # 4-region model
      ggfigs_ggpubr <- ggpubr::ggarrange( 
        # naive
        l_ggfigs[[ 1 ]][[ 1 ]], l_ggfigs[[ 2 ]][[ 1 ]], 
        l_ggfigs[[ 3 ]][[ 1 ]], l_ggfigs[[ 4 ]][[ 1 ]],
        # activated
        l_ggfigs[[ 1 ]][[ 2 ]], l_ggfigs[[ 2 ]][[ 2 ]], 
        l_ggfigs[[ 3 ]][[ 2 ]], l_ggfigs[[ 4 ]][[ 2 ]],
        # cd69+
        l_ggfigs[[ 1 ]][[ 3 ]], l_ggfigs[[ 2 ]][[ 3 ]], 
        l_ggfigs[[ 3 ]][[ 3 ]], l_ggfigs[[ 4 ]][[ 3 ]],
        labels = LETTERS[ 1 : ( 3 * length( pb.tissue ) ) ],
        ncol = length( pb.tissue ), nrow = 3 )  
    } else {
      # 3-region model
      ggfigs_ggpubr <- ggpubr::ggarrange( 
        # naive
        l_ggfigs[[ 1 ]][[ 1 ]], l_ggfigs[[ 2 ]][[ 1 ]], l_ggfigs[[ 3 ]][[ 1 ]], 
        # activated
        l_ggfigs[[ 1 ]][[ 2 ]], l_ggfigs[[ 2 ]][[ 2 ]], l_ggfigs[[ 3 ]][[ 2 ]], 
        # cd69+
        l_ggfigs[[ 1 ]][[ 3 ]], l_ggfigs[[ 2 ]][[ 3 ]], l_ggfigs[[ 3 ]][[ 3 ]], 
        labels = LETTERS[ 1 : ( 3 * length( pb.tissue ) ) ],
        ncol = length( pb.tissue ), nrow = 3 )  
    }
    
    ggsave( sprintf( "%s/%s.png", chain_dir, "all_trajectories" ), 
            ggfigs_ggpubr,
            width = 29.7 * 2.5, height = 21.0 * 2.5, units = "cm" )
    
    # 3-region model trajectories scaled to whole y axis
    if ( ( length( l_ggfigs_ext ) == 3 ) & ( length( pb.tissue ) == 3 ) )
    {
      ggfigs_ggpubr_ext <- ggpubr::ggarrange( 
        # naive
        l_ggfigs_ext[[ 1 ]][[ 1 ]], l_ggfigs_ext[[ 2 ]][[ 1 ]], 
        l_ggfigs_ext[[ 3 ]][[ 1 ]], 
        # activated
        l_ggfigs_ext[[ 1 ]][[ 2 ]], l_ggfigs_ext[[ 2 ]][[ 2 ]], 
        l_ggfigs_ext[[ 3 ]][[ 2 ]], 
        # cd69+
        l_ggfigs_ext[[ 1 ]][[ 3 ]], l_ggfigs_ext[[ 2 ]][[ 3 ]], 
        l_ggfigs_ext[[ 3 ]][[ 3 ]], 
        labels = LETTERS[ 1 : ( 3 * length( pb.tissue ) ) ],
        ncol = length( pb.tissue ), nrow = 3 )  
    ggsave( sprintf( "%s/%s.png", chain_dir, "all_trajectories_ext" ), 
            ggfigs_ggpubr_ext,
            width = 29.7 * 2.5, height = 21.0 * 2.5, units = "cm" )
    }
    
    if ( sel_type == "max" ) {
      # Trajectories are plotted for the best model fit
      model_id <- get_model_id( sel_models_file_name )
      fig_fname_start <- get_fig_fname_start( m = m,
        model_id = model_id, celltype = celltype, mcmc_pars = mcmc_pars )
      figs_traj_dir <- get_path( "figs_m_traj", m = m, ps = ps, celltype = celltype )
      ggsave( sprintf( "%s/%s_all_trajectories_%s.png",
                       figs_traj_dir, fig_fname_start, pb.tissue[ 2 ] ),
              ggfigs_ggpubr, 
              width = 29.7 * 2.5, height = 21.0 * 2.5, units = "cm" )
      
      gg_tissue_i_trajectories <- ggpubr::ggarrange( 
        l_ggfigs[[ 2 ]][[ 1 ]], 
        l_ggfigs[[ 2 ]][[ 2 ]], 
        l_ggfigs[[ 2 ]][[ 3 ]],
        ncol = 3, nrow = 1 )
      gg_tissue_i_trajectories %>% save( ., file = sprintf( 
        "%s/%s_%s_tis_trajs.rda",
        get_path( "model_dir", m = m ), fig_fname_start, pb.tissue[ 2 ] ) )
      
      if ( pb.tissue[ 2 ] == "spleen" ) {
        gg_tissue_blood_trajectories <- ggpubr::ggarrange( 
          l_ggfigs[[ 1 ]][[ 1 ]], 
          l_ggfigs[[ 1 ]][[ 2 ]], 
          l_ggfigs[[ 1 ]][[ 3 ]],
          ncol = 3, nrow = 1 )
        gg_tissue_blood_trajectories %>% save( ., file = sprintf( 
          "%s/%s_%s_tis_trajs.rda", 
          get_path( "model_dir", m = m ), fig_fname_start, "blood" ) )
      }
      
      # "Pooled others" trajectories:
      gg_tissue_i_trajectories <- ggpubr::ggarrange( 
        l_ggfigs[[ 3 ]][[ 1 ]], 
        l_ggfigs[[ 3 ]][[ 2 ]], 
        l_ggfigs[[ 3 ]][[ 3 ]],
        ncol = 3, nrow = 1 )
      gg_tissue_i_trajectories %>% save( ., file = sprintf( 
        "%s/%s_%s_pooled_trajs.rda",
        get_path( "model_dir", m = m ), fig_fname_start, pb.tissue[ 2 ] ) )
      
      # "Blood" trajectories:
      gg_tissue_i_trajectories <- ggpubr::ggarrange( 
        l_ggfigs[[ 1 ]][[ 1 ]], 
        l_ggfigs[[ 1 ]][[ 2 ]], 
        l_ggfigs[[ 1 ]][[ 3 ]],
        ncol = 3, nrow = 1 )
      gg_tissue_i_trajectories %>% save( ., file = sprintf( 
        "%s/%s_%s_blood_trajs.rda",
        get_path( "model_dir", m = m ), fig_fname_start, pb.tissue[ 2 ] ) )
      
      if ( ADD <- TRUE ) {
      for ( part_i in 1 : length( pb.tissue ) ) { 
        gg_tissue_i_trajectories <- ggpubr::ggarrange( 
          l_ggfigs[[ part_i ]][[ 1 ]], 
          l_ggfigs[[ part_i ]][[ 2 ]], 
          l_ggfigs[[ part_i ]][[ 3 ]],
          ncol = 3, nrow = 1 )
        gg_tissue_i_trajectories %>% save( ., file = sprintf( 
          "%s/%s_%s_part%i_trajs.rda", get_path( "model_dir", m = m ), 
          fig_fname_start, pb.tissue[ 2 ], part_i ) )
      }      
      
      }
    }
  }
}


test_eqeq <- function(
  m, parabio.fit = NA, 
  sel_chains, sel_type, iter = NA, write_csv = TRUE )
{
  model_name <- m$model_name
  celltype <- m$celltype; ps <- m$ps; 
  model_dir <- get_path( "model_dir", m = m )
  model_vars_dir <- get_path( "model_vars_dir", m = m )
  
  tol <- 1e-7
  chain_dir <- get_path( "chain_dir", m = m, model_dir = model_dir, 
                               sel_type = sel_type, sel_chains = sel_chains )
  
  if ( !is.na( iter ) ) {
    Q <- plot_par_densities_and_calc_Q(
      m = m,
      parabio.fit = parabio.fit, sel_type = sel_type, sel_chains = sel_chains, 
      iter = iter ) 
  } else {
    Q <- read_Q_matrix( chain_dir )  
  }    
  
  load( file = sprintf( "%s/other_vars.RData", model_vars_dir ) )
  load( file = sprintf( "%s/model_vars_.RData", model_vars_dir ) )
  
  eqeq_table <- tibble( 
    node = 1 : dim( Q )[ 1 ], influx = NA, outflux = NA, equal = NA, iter = iter )
  
  Qd0 <- Q - diag( diag( Q ) )
  for ( i in ( 1 : dim( Q )[ 1 ] ) ) {
    eqeq_table$outflux[ i ] <- -c( al_, al_ )[ i ] * Q[ i, i ]
    eqeq_table$influx[ i ] <- t( c( al_, al_ ) ) %*% Qd0[ , i ]
    eqeq_table$equal[ i ] <- abs( 
      eqeq_table$outflux[ i ] - eqeq_table$influx[ i ] ) < tol
  }
  if ( write_csv ) { write.csv( eqeq_table, sprintf( 
    "%s/aa_eqeq_table_iter_%s.csv", chain_dir,
    ifelse( is.na( iter), "all", iter ) ) ) }
  return( all( eqeq_table$equal ) )
}


get_Qcounts <- function(
    m, sel_chains, sel_type, 
    count_type = "Median" )
{
  celltype <- m$celltype; ps <- m$ps; 
  model_dir <- get_path( "model_dir", m = m )
  model_vars_dir <- get_path( "model_vars_dir", m = m )
  
  load( file = sprintf( "%s/pb.RData", model_vars_dir ) )
  # other_vars = parabio.data
  load( file = sprintf( "%s/other_vars.RData", model_vars_dir ) )
  # for al_
  load( file = sprintf( "%s/model_vars_.RData", model_vars_dir ) )
  chain_dir <- get_path( 
    "chain_dir", m = m, model_dir = model_dir, 
    sel_type = sel_type, sel_chains = sel_chains )
  
  d.celltype.counts <- read_csv( 
    sprintf( "%s/Total_counts/parabiosis_model_input_%s_counts.csv",
             ps$PROCESSED_PATH, celltype ) )
  f_al_to_ncells <- ( d.celltype.counts[ 
    d.celltype.counts$Tissue == "Blood", count_type ] %>%                       
      unlist %>% as.numeric ) / 
    sum( al_[ 1 : 3 ] )
  
  Q <- read_Q_matrix( chain_dir )
  
  Qcounts <- Q - diag( diag( Q ) )
  Qcounts_abs <- Qcounts
  for ( i in ( 1 : dim( Qcounts )[ 1 ] ) ) { 
    Qcounts_abs[ i, ] <- c( al_, al_ )[ i ] * Qcounts[ i, ] * f_al_to_ncells
    Qcounts[ i, ] <- c( al_, al_ )[ i ] * Qcounts[ i, ]
  }  
  diag( Qcounts ) <- NA
  diag( Qcounts_abs ) <- NA
  write.csv( Qcounts, file = sprintf( 
    "%s/Qcounts_matrix.csv", chain_dir ), quote = FALSE )
  write.csv( Qcounts_abs, file = sprintf(
    "%s/Qcounts_abs_matrix.csv", chain_dir ), quote = FALSE )
  
  return( Qcounts )
}


get_chain_output <- function( m, parabio.fit, 
                              sel_chains_table, i,
                              sel_models_file_name )
{
  Q <- plot_par_densities_and_calc_Q(
    m = m,
    parabio.fit = parabio.fit,
    sel_chains = sel_chains_table$sel_chains_set[ i ],
    sel_type = sel_chains_table$sel_chains_type[ i ] )
  
  
  final_fit_eqeq <- test_eqeq(
    m = m,
    sel_chains = sel_chains_table$sel_chains_set[ i ],
    sel_type = sel_chains_table$sel_chains_type[ i ] )
  
  
  # if debugging needed write where the equilibrium is not met
  if ( !final_fit_eqeq & sel_chains_table$sel_chains_set[ i ] == "1" & 
       sel_chains_table$sel_chains_type[ i ] == "each" ) {
    for ( iter in sample( 1 : 100, 10, replace = FALSE ) )
    {
      test_eqeq(
        m = m,
        parabio.fit = parabio.fit, 
        sel_chains = sel_chains_table$sel_chains_set[ i ],
        sel_type = sel_chains_table$sel_chains_type[ i ],
        iter = iter )           
    }
  }
  
  get_Qcounts(
    m = m,
    sel_chains = sel_chains_table$sel_chains_set[ i ],
    sel_type = sel_chains_table$sel_chains_type[ i ] )
  
  get_and_plot_trajectories(
    m = m,
    sel_chains = sel_chains_table$sel_chains_set[ i ],
    sel_type = sel_chains_table$sel_chains_type[ i ],
    sel_models_file_name = sel_models_file_name )  
}
