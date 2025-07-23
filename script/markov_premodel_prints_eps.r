# markov_premodel_prints_eps.r 

markov_premodel_prints_eps <- function( m = m, 
                                        vector_tissues,
                                        complem_tissues = c() )
{ 
  ps <- m$ps
  celltype <- m$celltype
  st <- m$st
  model_ver = st$MODEL_VER; 
  f_pb.week0i = st$PB_WEEK0I; f_eps = st$EPS;
  count_type <- "Median"
  
  figs_dir <- get_path( "figs_m", m = m )
  if ( !dir.exists( figs_dir ) ) { dir.create( figs_dir, recursive = TRUE ) }
  figs_traj_dir <- get_path( "figs_m_traj", m = m )
  if ( !dir.exists( figs_traj_dir ) ) { dir.create( figs_traj_dir, 
                                                    recursive = TRUE ) }
  figs_dwell_dir <- get_path( "figs_m_dwell", m = m )
  if ( !dir.exists( figs_dwell_dir ) ) { dir.create( figs_dwell_dir, 
                                                     recursive = TRUE ) }
  figs_rda_dir <- get_path( "figs_m_rda", m = m )
  if ( !dir.exists( figs_rda_dir ) ) { dir.create( figs_rda_dir, 
                                                   recursive = TRUE ) }
  figs_diag_dir <- get_path( "figs_m_diag", m = m )
  if ( !dir.exists( figs_diag_dir ) ) { dir.create( figs_diag_dir, 
                                                    recursive = TRUE ) }
  limit_dist_w12 = TRUE;
  
  get_normalized_proportion <- function( parabio.data.raw.i,
                                         cellstate, host_or_donor, p_norm )
  {
    parabio.data.raw.i_n <- parabio.data.raw.i
    host.celltype.cellstate <- sprintf( "%s.celltype.%s", "host", cellstate )
    donor.celltype.cellstate <- sprintf( "%s.celltype.%s", "donor", cellstate )
    
    # assign normalized values to column e.g. host.celltype.naive
    # if the sum of e.g. naive cells subset is > 0, i.e. DO NOT DIVIDE BY ZERO!       

    parabio.data.raw.i_n[ , host.celltype.cellstate ] <-
      ifelse( 
        ( ( parabio.data.raw.i[ , host.celltype.cellstate ] + 
              parabio.data.raw.i[ , donor.celltype.cellstate ] ) %>% 
            unlist %>% as.vector() ) > 0,   
        
        ( p_norm * 
            parabio.data.raw.i[ , host.celltype.cellstate ] / (
              parabio.data.raw.i[ , host.celltype.cellstate ] + 
                parabio.data.raw.i[ , donor.celltype.cellstate ] ) ) %>% 
          unlist() %>% as.numeric(), 
        
        p_norm * parabio.data.raw.i_n[ , "p_host" ]
        )
    
    # Now for donor:
    parabio.data.raw.i_n[ , donor.celltype.cellstate ] <-
      ifelse( 
        ( ( parabio.data.raw.i[ , host.celltype.cellstate ] + 
              parabio.data.raw.i[ , donor.celltype.cellstate ] ) %>% 
            unlist() %>% as.vector() ) > 0,   
        
        ( p_norm * parabio.data.raw.i[ , donor.celltype.cellstate ] / (
          parabio.data.raw.i[ , host.celltype.cellstate ] + 
            parabio.data.raw.i[ , donor.celltype.cellstate ] ) ) %>% 
          unlist %>% as.numeric(), 
        
        p_norm * ( 1 - parabio.data.raw.i_n[ , "p_host" ] )
        )
    return( parabio.data.raw.i_n )
  }
  
  normalize_tissue_states_in_time <- function( parabio.data.raw.i )
  {
    p_naive <- mean( unlist( parabio.data.raw.i[ , "host.celltype.naive" ] + 
                               parabio.data.raw.i[ , "donor.celltype.naive" ] ) )
    p_activ <- mean( unlist( parabio.data.raw.i[ , "host.celltype.activ" ] + 
                               parabio.data.raw.i[ , "donor.celltype.activ" ] ) )
    p_cd69p <- mean( unlist( parabio.data.raw.i[ , "host.celltype.cd69p" ] + 
                               parabio.data.raw.i[ , "donor.celltype.cd69p" ] ) )
    p_sum <- p_naive + p_activ + p_cd69p
    p_naive <- p_naive / p_sum
    p_activ <- p_activ / p_sum
    p_cd69p <- p_cd69p / p_sum
      
    parabio.data.raw.i$p_host <- ( parabio.data.raw.i[ , "host.celltype.naive" ] +
                                     parabio.data.raw.i[ , "host.celltype.activ" ] +
                                     parabio.data.raw.i[ , "host.celltype.cd69p" ] ) /
      ( parabio.data.raw.i[ , "host.celltype.naive" ] +
          parabio.data.raw.i[ , "host.celltype.activ" ] +
          parabio.data.raw.i[ , "host.celltype.cd69p" ] + 
          parabio.data.raw.i[ , "donor.celltype.naive" ] +
          parabio.data.raw.i[ , "donor.celltype.activ" ] +
          parabio.data.raw.i[ , "donor.celltype.cd69p" ] ) 
    parabio.data.raw.i_n <- 
      parabio.data.raw.i %>% 
      get_normalized_proportion( ., "naive", host_or_donor = NA, p_norm = p_naive ) %>% 
      get_normalized_proportion( ., "activ", host_or_donor = NA, p_norm = p_activ ) %>% 
      get_normalized_proportion( ., "cd69p", host_or_donor = NA, p_norm = p_cd69p )
      
    return( parabio.data.raw.i_n )
  }

  
    
  
  CORRECT_0 <- FALSE
  
  tissue.i1 <- vector_tissues[ 1 ]
  
  d.celltype.counts <- 
    read_csv( sprintf( "%s/Total_counts/parabiosis_model_input_%s_counts.csv",
                       ps$PROCESSED_PATH, celltype ) )
  
  
  if ( DEFINE_PARAMETERS_AND_SET_OPTIONS <- TRUE ) {
    
    data_dir <- get_path( "data", ps = ps, celltype = celltype ) 
    
    pb.week <- c( 0, 1, 2, 4, 8, 12 )
    if ( celltype == "B" ) {
      pb.week <- c( 0, 1, 2, 4, 8 )
    }
    
    l.data.file.name <- list()
    for( i in ( 1 : ( length( vector_tissues ) ) ) ) {
      l.data.file.name[[ i ]] <- sprintf( "%s_parabiosis_data_%s.csv", celltype,
                                          tolower( vector_tissues[ i ] ) ) 
    }
    if ( length( complem_tissues ) > 0 ) {
      l.complem.data.file.name <- list()
      for( i in ( 1 : ( length( complem_tissues ) ) ) ) {
        l.complem.data.file.name[[ i ]] <- sprintf( 
          "%s_parabiosis_data_%s.csv", celltype, tolower( complem_tissues[ i ] ) ) 
      }
    }
    
    vector_tissues_tmp <- vector_tissues
    if ( length( complem_tissues ) > 0 ) {
      vector_tissues_tmp <- c( vector_tissues, "Pooled__Others" )
    }
    
    pb.tissue <- c( "blood", tolower( vector_tissues_tmp ) )
    pb.tissue.label <- c( "Blood", vector_tissues_tmp )
    names( pb.tissue.label ) <- pb.tissue
    
    pb.hostdonor <- c( "host", "donor" )
    
    pb.parent.popul <- tolower( celltype )                                    
    pb.parent.popul.label <- celltype
    
    pb.sub.popul <- c( "naive", "activ", "cd69p" )
    pb.sub.popul.label <- c( "Naive", "Activated", "CD69+" )
    names( pb.sub.popul.label ) <- pb.sub.popul
    
    pb.popul <- paste( pb.parent.popul, pb.sub.popul, sep = "." )
    pb.popul.label <- paste( pb.parent.popul.label, pb.sub.popul.label, sep = " " )
    names( pb.popul.label ) <- pb.popul
    
    pb.week.n <- length( pb.week )
    pb.tissue.n <- length( pb.tissue )
    pb.popul.n <- length( pb.popul )

    if ( !is.na( f_pb.week0i ) ) { pb.week0i <- f_pb.week0i } else {
      pb.week0i <- 2
    }
    
    if ( !is.na( f_eps ) ) { 
      if ( pb.week0i != 0 ) {
      pb.week0i_eps <- f_pb.week0i - f_eps } else {
        pb.week0i_eps <- pb.week0i
        pb.week <- c( 0 + f_eps, pb.week[ -1 ] )
        CORRECT_0 <- TRUE
        }
    }
    
    pb.cell.ratio <- c( 750 )
    for ( i in ( 1 : ( length( vector_tissues ) ) ) ) {
      Blood_Count <- d.celltype.counts[ d.celltype.counts$Tissue == "Blood", count_type ]
      tissue.proportion.tissuei <- unlist( round( 
        d.celltype.counts[ 
          d.celltype.counts$Tissue == vector_tissues[ i ], count_type ] / 
          Blood_Count * 750, 2 ) )
      pb.cell.ratio <- c( pb.cell.ratio, tissue.proportion.tissuei )
    }
    if ( length( complem_tissues ) > 0 ) {                                      # adding artificial pool of tissues
      tissue.proportion.tissuei <- unlist( round(
        sum( d.celltype.counts[ 
          d.celltype.counts$Tissue %in% complem_tissues, count_type ] ) / 
        Blood_Count * 750, 2 ) )
      pb.cell.ratio <- c( pb.cell.ratio, tissue.proportion.tissuei )
    }
    
    names( pb.cell.ratio ) <- pb.tissue
    
    pb.figure.week.max <- 12

  }
  
  # set auxilary populations
  if ( SET_AUXILARY_POPULATIONS <- TRUE ) {
    pb.host.popul <- paste( rep( pb.hostdonor[ 1 ], each = length( pb.popul ) ),
                            pb.popul, sep = "." )
    
    pb.hostdonor.popul <- paste( rep( pb.hostdonor, each = length( pb.popul ) ),
                                 pb.popul, sep = "." )
    
    pb.hostdonor.tissue.popul <- paste(
      rep( pb.hostdonor, each = length( pb.tissue ) * pb.popul.n ),
      paste( rep( pb.tissue, each = pb.popul.n ), pb.popul, sep = "." ),
      sep = "." )
  }  
  
  
  # needed for the function estimate_week_alpha_simple_n()
  pb <- list(                                                                  
    pb.tissue, pb.tissue.label, pb.hostdonor, 
    pb.parent.popul, pb.parent.popul.label, 
    pb.sub.popul, pb.sub.popul.label, 
    pb.popul, pb.popul.label, 
    pb.week.n, pb.tissue.n, pb.popul.n,
    pb.week0i,
    pb.cell.ratio,
    pb.figure.week.max,
    pb.host.popul, pb.hostdonor.popul, pb.hostdonor.tissue.popul )
    names( pb ) <- 
      c( "pb.tissue", "pb.tissue.label", "pb.hostdonor",
         "pb.parent.popul", "pb.parent.popul.label",
         "pb.sub.popul", "pb.sub.popul.label",
         "pb.popul", "pb.popul.label",
         "pb.week.n", "pb.tissue.n", "pb.popul.n",
         "pb.week0i", "pb.cell.ratio", "pb.figure.week.max",
         "pb.host.popul", "pb.hostdonor.popul", "pb.hostdonor.tissue.popul" )
  
  # read and preprocess data
  if ( READ_AND_PREPROCESS_DATA <- TRUE ) {                                
  
    parabio.data.raw.1 <- read.csv( 
      file.path( data_dir, l.data.file.name[[ 1 ]] ) )
    names( parabio.data.raw.1 ) <- c( 
      "tissue", "week", "mouse", "pair",	
      "host.celltype.naive",	"host.celltype.activ",	"host.celltype.cd69p",	
      "donor.celltype.naive", "donor.celltype.activ",	"donor.celltype.cd69p" )
    if ( st$NORMALIZE )
    { 
      parabio.data.raw.blood <- parabio.data.raw.1[ 
        parabio.data.raw.1$tissue == "Blood", ]
      parabio.data.raw.blood_n <- normalize_tissue_states_in_time( 
        parabio.data.raw.i = parabio.data.raw.blood )
      
      parabio.data.raw.tissue.i1 <- parabio.data.raw.1[ 
        parabio.data.raw.1$tissue == tissue.i1, ]
      parabio.data.raw.tissue.i1_n <- normalize_tissue_states_in_time( 
        parabio.data.raw.i = parabio.data.raw.tissue.i1 )
      
      parabio.data.raw <- rbind( 
        parabio.data.raw.blood_n, parabio.data.raw.tissue.i1_n ) 
    } else {
      parabio.data.raw <- rbind( parabio.data.raw.1 )
    }
      
    # 4+ model parts
    if ( length( vector_tissues ) > 1 ) {
      for ( i in ( 2 : ( length( vector_tissues ) ) ) ) {
        parabio.data.raw.i <- read.csv( 
          file.path( data_dir, l.data.file.name[[ i ]] ) )
        names( parabio.data.raw.i ) <- c( 
          "tissue", "week", "mouse", "pair",	
          "host.celltype.naive",	"host.celltype.activ",	"host.celltype.cd69p",	
          "donor.celltype.naive", "donor.celltype.activ",	"donor.celltype.cd69p" )
        
        parabio.data.raw.i <- parabio.data.raw.i[ parabio.data.raw.i$tissue != "Blood", ]
        if ( st$NORMALIZE ) { 
          parabio.data.raw.i_n <- 
            normalize_tissue_states_in_time( parabio.data.raw.i )
          parabio.data.raw <- rbind( parabio.data.raw, parabio.data.raw.i_n )
        } else {
        parabio.data.raw <- rbind( parabio.data.raw, parabio.data.raw.i )
        }
      }
    }
    
    Blood_Count <- d.celltype.counts[ d.celltype.counts$Tissue == "Blood", count_type ]
    
    if ( length( complem_tissues ) > 0 ) {
      parabio.data.raw.complem <- NULL
      for ( j in ( 1 : ( length( complem_tissues ) ) ) ) { 
        parabio.data.raw.j <- read.csv( file.path( data_dir, l.complem.data.file.name[[ j ]] ) )
        names( parabio.data.raw.j ) <- c( 
          "tissue", "week", "mouse", "pair",	
          "host.celltype.naive",	"host.celltype.activ",	"host.celltype.cd69p",	
          "donor.celltype.naive", "donor.celltype.activ",	"donor.celltype.cd69p" )
        parabio.data.raw.j <- parabio.data.raw.j[ parabio.data.raw.j$tissue != "Blood", ]
        parabio.data.raw.j$total.count <- 
          d.celltype.counts[ d.celltype.counts$Tissue == complem_tissues[ j ], 
                             count_type ] %>% unlist
          
        parabio.data.raw.complem <- rbind( parabio.data.raw.complem, 
                                           parabio.data.raw.j )
      }
      parabio.data.raw.complem.pooled <- 
        parabio.data.raw.complem %>% as_tibble() %>% 
        group_by( mouse ) %>%
        summarise( 
          tissue = "Pooled__Others",
          week = as.integer( mean( week ) ),
          pair = as.integer( mean( pair ) ),
          host.celltype.naive = weighted.mean( x = host.celltype.naive, w = total.count ),
          host.celltype.activ = weighted.mean( x = host.celltype.activ, w = total.count ),
          host.celltype.cd69p = weighted.mean( x = host.celltype.cd69p, w = total.count ),
          donor.celltype.naive = weighted.mean( x = donor.celltype.naive, w = total.count ),
          donor.celltype.activ = weighted.mean( x = donor.celltype.activ, w = total.count ),
          donor.celltype.cd69p = weighted.mean( x = donor.celltype.cd69p, w = total.count ) )
        
      if ( st$NORMALIZE ) { 
        parabio.data.raw.complem.pooled_n <-
          normalize_tissue_states_in_time( parabio.data.raw.complem.pooled )
        parabio.data.raw <- rbind( parabio.data.raw, parabio.data.raw.complem.pooled_n )
      } else {
        parabio.data.raw <- rbind( parabio.data.raw, parabio.data.raw.complem.pooled )
      }
    }
    
    
    parabio.data.raw <- parabio.data.raw[ rowSums( parabio.data.raw[ , 5:10 ] ) != 0, ]  

    parabio.data <- parabio.data.raw
    
    parabio.data$tissue <- tolower( parabio.data$tissue )
    
    stopifnot( parabio.data$tissue %in% pb.tissue )
    
    if ( CORRECT_0 ) {
      parabio.data$week[ parabio.data$week == 0 ] <- 0 + f_eps 
    }
    stopifnot( parabio.data$week %in% pb.week )
  }
  
  
  if ( FIT_DIRICHLET_AT_W0_AND_GET_DATA_FOR_LATER <- TRUE ) {
    if ( !CORRECT_0 ) {
    week00.alpha <- estimate_week_alpha_simple_n( 
      parabio.data, w = 0, n.tis = pb.tissue.n, pb = pb )
    week0i.alpha <- estimate_week_alpha_simple_n( 
      parabio.data, w = pb.week0i, n.tis = pb.tissue.n, pb = pb ) 
    } else {
        week00.alpha <- estimate_week_alpha_simple_n( 
          parabio.data, w = 0 + f_eps, n.tis = pb.tissue.n, pb = pb )
        week0i.alpha <- estimate_week_alpha_simple_n( 
          parabio.data, w = pb.week0i + f_eps, n.tis = pb.tissue.n, pb = pb )
      }

    
    # get data for week > initial week
    parabio.model.week <- parabio.data[ parabio.data$week > pb.week0i_eps, "week" ]
    
    parabio.model.tissue <- parabio.data[ parabio.data$week > pb.week0i_eps, "tissue" ]
    
    parabio.model.x <- data.matrix( parabio.data[ 
      parabio.data$week > pb.week0i_eps,
      gsub( 
        pattern = "\\.(.*)\\.", 
        replacement = ".celltype.", 
        pb.hostdonor.popul ) ] )
    stopifnot( parabio.model.x >= 0 )
    parabio.model.x <- sweep( parabio.model.x, 1, rowSums( parabio.model.x ), "/" )
    parabio.model.x <- model_x_remove_zeros( parabio.model.x )
    
    parabio.model.x <- sweep( parabio.model.x, 1, rowSums( parabio.model.x ), "/" )
    parabio.model.x 
  }
  
  
  # set variables for markov model
  if ( SET_VARIABLES_FOR_MARKOV_MODEL <- TRUE ) {
    if ( !is.na( st$HOST_DONOR_RATE_ ) ) { 
      host_donor_rate_ <- st$HOST_DONOR_RATE_ %>% as.numeric } else { 
        host_donor_rate_ <- 100 }
    if ( !is.na( st$USE_HDR2_ ) ) {
      use_hdr2_ <- st$USE_HDR2_ %>% as.integer() } else {
        use_hdr2_ <- 0L
      }
    if ( !is.na( st$BRF ) ) {
      brf <- st$BRF %>% as.numeric() } else {
        brf <- 1
      }
    if ( !is.na( st$MAX_FLOW_TO_N3_ ) ) {
      max_flow_to_N3_ <- st$MAX_FLOW_TO_N3_ %>% as.numeric() } else {
        max_flow_to_N3_ <- 1    
      }
    
    n_ <- nrow( parabio.model.x )
    wn_ <- sum( pb.week > pb.week0i_eps )
    tn_ <- pb.tissue.n
    pn_ <- pb.popul.n
    
    w0i_ <- pb.week0i_eps
    w_ <- pb.week[ pb.week > pb.week0i_eps ]
    
    wi_ <- match( parabio.model.week, w_ )
    ti_ <- match( parabio.model.tissue, pb.tissue )
    x_ <- parabio.model.x
    
    al_ <- as.vector( sapply( pb.tissue, function( tis )
      ( pb.cell.ratio[ tis ] / sum( pb.cell.ratio ) ) *
        ( week00.alpha[ tis, ] / sum( week00.alpha[ tis, ] ) )
    ) )
    
    al_old_ <- al_
    if ( limit_dist_w12 ) {
      if ( celltype != "B" ) { w12 <- 12 } else { w12 <- 8 }
      week12.alpha <- estimate_week_alpha_simple_n( 
        parabio.data, w = w12, n.tis = pb.tissue.n, pb = pb )
      
      al12_tmp_ <- as.vector( sapply( pb.hostdonor, function( hd )
        sapply( pb.tissue, function( tis )
          ( pb.cell.ratio[ tis ] / sum( pb.cell.ratio ) ) *
            ( week12.alpha[ tis,
                            grep( sprintf( "^%s\\.", hd ), colnames( week12.alpha ) ) ] /
                sum( week12.alpha[ tis, ] ) )
        )
      ) )
      
      al12_ <- al12_tmp_[ 1 : ( pn_ * tn_ ) ] + 
        al12_tmp_[ ( 1 + pn_ * tn_ ) : ( 2 * pn_ * tn_ ) ]
      al_ <- al12_
    }
    
    {
      if ( pb.week0i == 0 )
        a0i_ <- c( al_, rep( 0, pb.tissue.n * pb.popul.n ) )
      else
        a0i_ <- as.vector( sapply( pb.hostdonor, function( hd )
          sapply( pb.tissue, function( tis )
            ( pb.cell.ratio[ tis ] / sum( pb.cell.ratio ) ) *
              ( week0i.alpha[ tis,
                              grep( sprintf( "^%s\\.", hd ), colnames( week0i.alpha ) ) ] /
                  sum( week0i.alpha[ tis, ] ) )
          )
        ) )
    }
    
    w_a <- f_pb.week0i
    if ( CORRECT_0 ) { w_a <- f_pb.week0i + f_eps }
    w_b <- min( ( pb.week )[ pb.week > f_pb.week0i ] )  # the next after start
    if ( NEW_2_4 <- TRUE ) {
      # week a
      week0a.alpha <- estimate_week_alpha_simple_n( 
        parabio.data, w = w_a, n.tis = pb.tissue.n, pb = pb )
      
      # week b
      week0b.alpha <- estimate_week_alpha_simple_n( 
        parabio.data, w = w_b, n.tis = pb.tissue.n, pb = pb )
      
      a04_ <- as.vector( sapply( pb.hostdonor, function( hd )
        sapply( pb.tissue, function( tis )
          ( pb.cell.ratio[ tis ] / sum( pb.cell.ratio ) ) *
            ( week0b.alpha[ tis,
                            grep( sprintf( "^%s\\.", hd ), colnames( week0b.alpha ) ) ] /
                sum( week0b.alpha[ tis, ] ) )
        )
      ) )
      # extrapolation:
      if ( !is.na( f_eps ) ) {
      week0a.alpha_eps_tmp <- week0a.alpha - ( week0b.alpha - week0a.alpha ) * f_eps
      week0a.alpha_eps_tmp[ week0a.alpha_eps_tmp <= 0 ] <- 
        week0a.alpha[ week0a.alpha_eps_tmp <= 0 ]
      a02_eps_ <- as.vector( sapply( pb.hostdonor, function( hd )
        sapply( pb.tissue, function( tis )
          ( pb.cell.ratio[ tis ] / sum( pb.cell.ratio ) ) *
            ( week0a.alpha_eps_tmp[ tis,
                            grep( sprintf( "^%s\\.", hd ), colnames( week0a.alpha_eps_tmp ) ) ] /
                sum( week0a.alpha_eps_tmp[ tis, ] ) )
        )
      ) )
      }
    }
    if ( !is.na( f_eps ) ) {
      a0i_ <- a02_eps_
    }
  }
    
  names( a0i_ ) <- pb.hostdonor.tissue.popul
    
  # save model pb and starting.data_ to files and .RData - to create a function
  if ( SAVING_FOR_LOG.LIK_MANUAL_CALC <- TRUE ) {
    model_vars_dir <- get_path( 
      type = "model_vars_dir", m = m, 
      celltype = celltype, tissue = tissue.i1, ps = ps, mcmc_pars = m$mcmc_pars, 
      model_name = m$model_name, model_ver = model_ver )
  
    if ( !dir.exists( model_vars_dir ) ) 
      dir.create( model_vars_dir, recursive = TRUE )
    
    save( pb.tissue, pb.tissue.label, pb.hostdonor, 
          pb.parent.popul, pb.parent.popul.label, 
          pb.sub.popul, pb.sub.popul.label, 
          pb.popul, pb.popul.label, 
          pb.week, pb.week.n, pb.tissue.n, pb.popul.n,
          pb.week0i, pb.cell.ratio, pb.figure.week.max,
          pb.host.popul, pb.hostdonor.popul, pb.hostdonor.tissue.popul,
          file = sprintf( "%s/pb.RData", model_vars_dir ) )
    write.csv( pb.cell.ratio, row.names = TRUE,
               sprintf( "%s/pb_cell_ratio.csv", model_vars_dir ) )

    # e.g., n_cells_born_per_day is 100000 for Treg
    d.cells.born.counts <- 
      read_csv( sprintf( "%s/Total_counts/cells_born_counts.csv",
                         ps$PROCESSED_PATH ) )
    n_born_cells_day_celltype <- d.cells.born.counts[ 
      d.cells.born.counts$celltype == celltype, "n_cells_born_per_day" ] %>% 
      unlist %>% as.numeric

    d.celltype.counts <- 
      read_csv( sprintf( "%s/Total_counts/parabiosis_model_input_%s_counts.csv",
                         ps$PROCESSED_PATH, celltype ) )
    n_blood_cells <- d.celltype.counts[ 
      d.celltype.counts$Tissue == "Blood", count_type ] %>% 
      unlist %>% as.numeric
    b_rate_rel_to_blood_naive_ <- brf * 
      n_born_cells_day_celltype * sum( al_[ 1:3 ] ) / ( n_blood_cells * al_[ 1 ] )

    
save( host_donor_rate_, n_, wn_, tn_, pn_, w0i_, w_, wi_, ti_, x_, al_, a0i_, 
      b_rate_rel_to_blood_naive_, 
      # use_wghs_, 
      use_hdr2_,
      max_flow_to_N3_, 
      file = sprintf( "%s/model_vars_.RData", model_vars_dir ) )

    save( parabio.data,
          file = sprintf( "%s/other_vars.RData", model_vars_dir ) )
    write.csv( a0i_, sprintf( "%s/a0i_.csv", model_vars_dir ),
               row.names = names( a0i_ ) )
    write.csv( x_, row.names = TRUE,
               sprintf( "%s/x_.csv", model_vars_dir ) )
    write.csv( data.frame( wi_, parabio.model.week ), row.names = TRUE,
               sprintf( "%s/wi_.csv", model_vars_dir ) )
    write.csv( ti_, row.names = TRUE,
               sprintf( "%s/ti_.csv", model_vars_dir ) )
    write.csv( al_, row.names = TRUE,
               sprintf( "%s/al_.csv", model_vars_dir ) )
    write.csv( al12_tmp_, row.names = TRUE,
               sprintf( "%s/al12_tmp_.csv", model_vars_dir ) )
    write.csv( al_old_, row.names = TRUE,
               sprintf( "%s/al_old_.csv", model_vars_dir ) )
  }
}

