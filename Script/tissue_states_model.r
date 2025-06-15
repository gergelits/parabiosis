# tissue_states_model.r

get_tissue_states_model_exit <- function( Q, al_ )
{
  p_i <- al_[ 4 : 6 ] / sum( al_[ 4 : 6 ] )
  q45 <- Q[ 4, 5 ]; q41 <- Q[ 4, 1 ]
  q56 <- Q[ 5, 6 ]; q52 <- Q[ 5, 2 ]; q51 <- Q[ 5, 1 ]
  q65 <- Q[ 6, 5 ]; q63 <- Q[ 6, 3 ]; q61 <- Q[ 6, 1 ]
  q4 <- q45 + q41
  q5 <- q56 + q52 + q51
  q6 <- q65 + q61 + q63
  p4_in <- q45 / q4
  p5_in <- q56 / q5
  p6_in <- q65 / q6
  
  p4_die <- 0 # q41 / q4 / 2
  p5_die <- q51 / q5
  p6_die <- q61 / q6

  p4_bl <- q41 / q4
  p5_bl <- q52 / q5
  p6_bl <- q63 / q6
  
  m456 <- tibble( p_i = p_i,
                  state = 4 : 6,
                  exp_rate = c( q4, q5, q6 ),
                  p_in = c( p4_in, p5_in, p6_in ),
                  p_die = c( p4_die, p5_die, p6_die ),
                  p_bl = c( p4_bl, p5_bl, p6_bl ),
                  if_stay_goto = c( 5, 6, 5 ) )
  return( m456 )
}


get_tissue_states_model_entry <- function( Q, al_ )
{
  q14 <- Q[ 1, 4 ]
  q25 <- Q[ 2, 5 ]; q45 <- Q[ 4, 5 ]; q65 <- Q[ 6, 5 ]
  q36 <- Q[ 3, 6 ]; q56 <- Q[ 5, 6 ]
  q41 <- Q[ 4, 1 ]; q52 <- Q[ 5, 2 ]; q51 <- Q[ 5, 1 ]; 
  q63 <- Q[ 6, 3 ]; q61 <- Q[ 6, 1 ]

  q4 <- q45 + q41
  q5 <- q56 + q52 + q51
  q6 <- q65 + q61 + q63

  s4_entry <- al_[ 1 ] * q14 
  s5_entry <- al_[ 2 ] * q25 + al_[ 4 ] * q45 + al_[ 6 ] * q65
  s6_entry <- al_[ 3 ] * q36 + al_[ 5 ] * q56
  
  p4_out <- al_[ 1 ] * q14 / s4_entry
  p5_out <- al_[ 2 ] * q25 / s5_entry
  p6_out <- al_[ 3 ] * q36 / s6_entry
  
  p5_in_s4 <- al_[ 4 ] * q45 / s5_entry
  p5_in_s6 <- al_[ 6 ] * q65 / s5_entry
  p5_in <- p5_in_s4 + p5_in_s6
  p6_in <- al_[ 5 ] * q56 / s6_entry
  
  come_from_s <- list()
  come_from_s[[ "s4" ]] <- c( "out" )
  come_from_s[[ "s5" ]] <- c( "out", "s4", "s6" )
  come_from_s[[ "s6" ]] <- c( "out", "s5" )
  come_from_p <- list()
  come_from_p[[ "s4" ]] <- c( p4_out )
  come_from_p[[ "s5" ]] <- c( p5_out, p5_in_s4, p5_in_s6 )
  come_from_p[[ "s6" ]] <- c( p6_out, p6_in )
  exp_rate <- c( q4, q5, q6 ); names( exp_rate ) <- c( "s4", "s5", "s6" )
  
  entry <- list( come_from_s, come_from_p, exp_rate )
  names( entry ) <- c( "come_from_s", "come_from_p", "exp_rate" )
  return( entry )
}


sim_tissue_states_dwell_time_exit <- function( n_sim = 1000, m456, s_init = 0 )
{
  stopifnot( s_init %in% c( 0, 4, 5, 6 ) )

  sim_steps <- tibble( i_sim = NULL, state = NULL, occur = NULL,
                       t = NULL, time = NULL )
  
  for ( i in 1 : n_sim ) {
    if ( s_init == 0 ) {
      s_curr <- sample( x = m456$state, size = 1, prob = m456$p_i )
    } else {
      s_curr <- s_init
    }
    # coding of "step": 0 = naive; odd = activ; even = cd69p
    step <- s_curr - 4 
    t <- 0
    
    cell_next <- "in"
    while ( cell_next == "in" ) {
      curr_dwell_time <- 
        rexp( n = 1, rate = m456$exp_rate[ m456$state == s_curr ] )
      
      sim_step_curr <- tibble( 
        i_sim = i, 
        state = sprintf( "s%s", s_curr ),
        occur = sprintf( "o%s", stringr::str_pad( step, 3, pad = "0" ) ),
        t = t,
        time = curr_dwell_time )
      sim_steps <- bind_rows( sim_steps, sim_step_curr )
      
      cell_next <- sample( x = c( "in", "die", "blood" ), size = 1,
                           prob = c( m456$p_in[ m456$state == s_curr ],
                                     m456$p_die[ m456$state == s_curr ],
                                     m456$p_bl[ m456$state == s_curr ] ) )
  
      if ( cell_next == "in" ) {
        s_curr <- m456$if_stay_goto[ m456$state == s_curr ]
        step <- step + 1
        t <- t + 1
      } else {
        exit_type <- ifelse( 
          s_curr != 4, sprintf( "%s_s%s", cell_next, s_curr ), "exit_s4" )
        
        sim_exit_curr <- tibble( i_sim = i, 
                                 state = exit_type,
                                 occur = "exit",
                                 t = t + 1,
                                 time = 0 )
        sim_steps <- bind_rows( sim_steps, sim_exit_curr )
      }
    }
  }
  
  return( sim_steps )
}


sim_tissue_states_past_future <- function( m456_en, sim_future )
{
  d_init_states <- sim_future %>% 
    group_by( i_sim ) %>% 
    filter( row_number() == 1 ) %>% 
    dplyr::select( i_sim, state )
  sim_steps_p <- tibble( i_sim = NULL, state = NULL, occur = NULL, 
                         t = NULL, time = NULL )
  
  for ( i in d_init_states$i_sim )
  {
    ss_curr <- d_init_states$state[ i ]; 
    cell_prev <- sample( x = m456_en$come_from_s[[ ss_curr ]], size = 1, 
                         prob = m456_en$come_from_p[[ ss_curr ]] )
    t <- -1
    while ( cell_prev != "out" ) {
      prev_dwell_time <- rexp( n = 1, rate = m456_en$exp_rate[ cell_prev ] )
      sim_step_prev <- tibble( i_sim = i, 
                               state = cell_prev,
                               occur = NA,
                               t = t,
                               time = prev_dwell_time )
      sim_steps_p <- bind_rows( sim_steps_p, sim_step_prev )
      
      cell_prev <- sample( x = m456_en$come_from_s[[ cell_prev ]], size = 1, 
                           prob = m456_en$come_from_p[[ cell_prev ]] )
      t <- t - 1
    }
  }

  sim_past_fut <-
    sim_future %>% 
    mutate( from = "future" ) %>% 
    bind_rows( ., sim_steps_p %>% mutate( from = "past" ) ) %>% 
    arrange( i_sim, t ) %>% 
    group_by( i_sim ) %>%
    mutate( t_pos = t - min( t ) ) %>% 
    mutate( state_start = ifelse( t == min( t ), substring( state, 2 ), -1 ) ) %>%
    mutate( state_start = max( as.numeric( state_start ) ) ) %>% 
    mutate( occur_old = occur ) %>% 
    mutate( occur = ifelse( 
      ! is.na( occur_old ) & occur_old == "exit", "exit", sprintf( 
        "o%s", stringr::str_pad( state_start - 4 + t_pos, 3, pad = "0" ) ) ) ) %>% 
    ungroup()
  
  return( sim_past_fut )
}




