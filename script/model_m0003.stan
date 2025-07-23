// 18 states
data {
  real host_donor_rate_;   // hdr can be set as fixed - used for testing only
  
  int n_;           // number of parabiosis samples
  
  int wn_;          // number of weeks greater than initial week
  int tn_;          // number of tissues/model compartments (3)                                             
  int pn_;          // number of cell populations (3)
  
  real w0i_;        // initial week
  array[ wn_] real w_;   // actual values of weeks greater than the initial week
  
  int wi_[ n_ ];    // index of sample weeks
  int ti_[ n_ ];    // index of sample tissues
  array[ n_ ] simplex[ 2*pn_ ] x_;    // samples
  
  // limiting distribution of single mouse
  simplex[ tn_*pn_ ] al_;
  
  // distribution of initial week
  simplex[ 2*tn_*pn_ ] a0i_;

  // passed model parameters
  real b_rate_rel_to_blood_naive_;
  int use_hdr2_;
  real max_flow_to_N3_; 
}


transformed data {
  int d = 2*tn_*pn_;                   // number of populations per paired mice
  // upper bounds
  real q65_u = 10;
  real q63_u = 0.10 / max( { al_[6], al_[9] } );
}


parameters {
  // free model parameters
  // generally, qij represents flow rate from model node i (Ni) to model node j (Nj)
  // model nodes: N1 = blood naive,  N2 = blood activated,  N3 = blood CD69+
  //              N4 = tissue naive, N5 = tissue activated, N6 = tissue CD69+
  // arrays/vectors such as q63_i represent: 
  //    q63_i[1] = flow from N6[1] to N3, where N6[1] is N6 above
  //    q63_i[2] = flow from N6[2] to N3, where N6[2] is CD69+ in combined pool of other tissue

  array[ tn_-1 ] real< lower=0, upper=q63_u > q63_i;
  
  // free helper parameters - allows to uniformly sample within the range
  // predefined by precalculated constrains (to also fulfill total_entry = total_exit)
  array[ tn_-1 ] real< lower=0, upper=1 > q14_base_i;
  array[ tn_-1 ] real< lower=0, upper=1 > q36_base_i;
  array[ tn_-1 ] real< lower=0, upper=1 > q56_base_i;
  array[ tn_-1 ] real< lower=0, upper=1 > q65_base_i;
  array[ tn_-1 ] real< lower=0, upper=1 > q52_base_i;
  array[ tn_-1 ] real< lower=0, upper=1 > q25_base_i;
  array[ tn_-1 ] real< lower=0, upper=1 > q51_base_i;
  
  // small CD69+ de-differentation in blood
  real< lower=0, upper=0.010 > q32;                
  
  // birth rate factor
  real< lower=1, upper=30 > brf_free_param;        
  
  // 3 free between mouse parameters
  real< lower=3, upper=1000 > host_donor_rate;     // naive host-donor exchange
  real< lower=0, upper=200 > host_donor_rate2;     // activated host-donor exchange
  real< lower=0, upper=3000 > host_donor_rate3; // CD69+ host-donor exchange
  
  // 15 nuisance parameters 
  // 5 x 3 (5 weeks: 1, 2, 4, 8, 12) x (3 tissues: tissue_i, blood, others)
  // minus 1 x 3 (week 1 parameters set as fixed to request fit through blood)
  array[ wn_, tn_ ] real< lower=0, upper=1000 > af0;  
}

//// DO NOT DELETE the following comments about definitions in R

// start of definitions for R   
// af = af0 = matrix( nrow = wn_, ncol = tn_ )
// q14_min_I_i = q14_max_I_i = q14_base_i = numeric( tn_-1 )
// q41_u_code_i = q45_u_code_i = numeric( tn_-1 )
// al1_q14_min_i = al4_q41_flow_i = numeric( tn_-1 )
// q51_u_code_i = q61_u_code_i = numeric( tn_-1 )
// q14_simplex_i = numeric( tn_-1 )
// wgh4_free_i = wgh5_free_i = wgh6_free_i = numeric( tn_-1 )
// q41_flow_i = numeric( tn_-1 )
// wgh4_i = wgh5_i = wgh6_i = numeric( tn_-1 )
// q36_propos_i = numeric( tn_-1 )
// q25_propos_i = numeric( tn_-1 )
// q14_propos_i = numeric( tn_-1 )
// q41_base_i = numeric( tn_-1 )
// q14_i = q14_min_i = numeric( tn_-1 )
// q45_max_i = numeric( tn_-1 )
// q41_i = q51_i = q52_i = q61_i = q63_i = q65_i = q56_base_i = numeric( tn_-1 ) 
// q41a_i = q41b_i = q56_i = q45_i = numeric( tn_-1 )                          
// q44_i = q55_i = q66_i = q36_i = q25_i = numeric( tn_-1 )
// q36_simplex_i = q36_sum_max_simplex_i = q36_sum_min_simplex_i = numeric( tn_-1 )
// q25_simplex_i = q25_sum_max_simplex_i = q25_sum_min_simplex_i = numeric( tn_-1 )
// q23_min = q23_max = numeric( tn_-1 )
// q61_q65_sum_i = q66_min_i = q66_max_i = q65_base_i = numeric( tn_-1 )
// q61_q65_sum_propos_i = numeric( tn_-1 )
// q36_max_S_i = q36_max_I_i = q36_max_i = q36_min_S_i = q36_min_I_i = q36_min_i = q36_base_i = numeric( tn_-1 )
// q25_max_S_i = q25_max_I_i = q25_max_i = q25_min_S_i = q25_min_I_i = q25_min_i = q25_base_i = numeric( tn_-1 )
// q25_inminmax_simplex_i = numeric( tn_-1 )
// q25_S_i = numeric( tn_-1 )
// al_2_q25__al_1 = numeric( tn_-1 )
// q51_min_i = q51_max_i = q51_base_i = numeric( tn_-1 )
// q25_min_N5_i = q25_max_N5_i = q52_min_N5_i = q52_max_N5_i = q52_min_N2_i = q52_max_N2_i = q52_min_i = q52_max_i = q52_base_i = numeric( tn_-1 )

// end of definitions for R

transformed parameters {
  array[ wn_, tn_ ] real< lower=0, upper=1000 > af;  
  
  array[ tn_-1 ] real< lower=0 > q14_propos_i; // for testing to avoid overflow
  array[ tn_-1 ] real< lower=0 > q14_max_I_i;
  
  array[ tn_-1 ] real< lower=0 > q14_i;
  array[ tn_-1 ] real< lower=0 > q14_min_I_i;
  
  array[ tn_-1 ] real< lower=0 > q56_i;
  array[ tn_-1 ] real< lower=0 > q36_min_I_i;
  array[ tn_-1 ] real< lower=0 > q36_max_I_i;
  array[ tn_-1 ] real< lower=0 > q36_i;

  real< lower=0 > q23;

  array[ tn_-1 ] real< lower=0 > q51_min_i;
  array[ tn_-1 ] real< lower=0 > q51_max_i;
  array[ tn_-1 ] real< lower=0 > q51_i;
  array[ tn_-1 ] real< lower=0 > q45_max_i;
  array[ tn_-1 ] real< lower=0 > q45_i;
  array[ tn_-1 ] real< upper=0 > q44_i;
  array[ tn_-1 ] real< upper=0 > q55_i;

  array[ tn_-1 ] real< upper=0 > q66_min_i;
  array[ tn_-1 ] real< upper=0 > q66_max_i;
  array[ tn_-1 ] real< upper=0 > q66_i;
  array[ tn_-1 ] real< lower=0 > q61_q65_sum_propos_i;
  array[ tn_-1 ] real< lower=0 > q61_q65_sum_i;
  array[ tn_-1 ] real< lower=0 > q61_i;
  array[ tn_-1 ] real< lower=0, upper=q65_u > q65_i;          
  array[ tn_-1 ] real< lower=0 > q41_i;
  real< lower=0 > q12_tmp;
  
  real q17 = host_donor_rate;  
  real q28 = host_donor_rate;
  real q39 = host_donor_rate;
  
  // Inconveniently, the following 7 values need to be set twice. Here and 
  // below after definition of "nonfree" variables
  real q12_l_code = -1000;
  real q12_u_code = -1000;
  real q41_u_code = -1000;
  real q52_u_code = -1000;
  real q56_q65_u_code = -1000;
  real q23_u_code = -1000;
  real q1j_min = -1000;
  
  real< lower=0 > q12;
  real< upper=0 > q33;
  real< upper=0 > q22;
  real< upper=0 > q11;
  
  //// real sum_al_4_q41_flow;
  real sum_q1j;
  
  real< lower=0 > q14_sum_max;
  real< lower=0 > q14_sum_min;
  real< lower=0 > q14_min_last;
  
  real sum_al_6_q63;
  real sum_al_5_q52;
  real sum_al_5_q52_max_N5;
  real sum_al_5_q52_min_N5;
  
  
  array[ tn_-1 ] real< lower=0 > q36_propos_i;
  real< lower=0 > q36_min_last;
  real< lower=0 > q36_max_last;
  
  // constrains to be calculated for the q25 and q52 parameters
  array[ tn_-1 ] real< lower=0 > q25_min_N5_i;
  array[ tn_-1 ] real< lower=0 > q25_max_N5_i;
  array[ tn_-1 ] real< lower=0 > q25_i;
  
  array[ tn_-1 ] real< lower=0 > q52_min_N5_i;
  array[ tn_-1 ] real< lower=0 > q52_max_N5_i;
  array[ tn_-1 ] real< lower=0 > q52_min_N2_i;
  array[ tn_-1 ] real< lower=0 > q52_max_N2_i;
  array[ tn_-1 ] real< lower=0 > q52_min_i;
  array[ tn_-1 ] real< lower=0 > q52_max_i;
  array[ tn_-1 ] real< lower=0 > q52_i;
  real< lower=0 > al2_out;
  
  array[ tn_-1 ] real q41_u_code_i;
  array[ tn_-1 ] real q45_u_code_i;
  
  array[ tn_-1 ] real q61_u_code_i;
  array[ tn_-1 ] real q51_u_code_i;
  
  
  // !!! DO NOT REWRITE THE FOLLOWING LINE - IT IS DETECTED BY R PARSER
  // start of nonfree variables
    q12_l_code = 0;
    q12_u_code = 0.01;
    q41_u_code = 3;
    q52_u_code = 10;
    q56_q65_u_code = 10;
    q23_u_code = 0.010;
    q1j_min = 0.0000001;

      
    {{//REMOVE IN R                                                             //// The curly brackets "{{" are allowing me to have local variables in "trnansformed parameters {...}" section
      //// !! WARNING: To comply with my current Stan->R parser, the following lines cannot by wrapped 
      //// TODO: delete the next line - MORE COMPLICATED
      int n_t_wgh_1 = 1;                                                        //// number of tissues / model parts without death weighting (i.e. wgh[ 1:n_t_wgd_1 ] = 1)
      int of = 3 * tn_-1; int ho = 3 * tn_;                                     //// hack for the R parser
      real q36_sum_max = -1000;                                                 //// Defined here to get parsed into R environment
      real q36_sum_min = -1000;
      
      //// It is possible to predefine a fixed hdr value:
      if ( host_donor_rate_ > 0 ) {{ q17 = host_donor_rate_; q28 = q17; q39 = q17;    }}
      //// Should separate parameters for host_donor_rate2, host_donor_rate3 be used?
      if ( use_hdr2_ == 1 )       {{ q28 = host_donor_rate2; q39 = host_donor_rate3;  }}

      for ( i in 1                 : n_t_wgh_1 ) {{ q41_u_code_i[ i ] = 0.01;  q45_u_code_i[ i ] = 5;   }}
      for ( i in ( n_t_wgh_1 + 1 ) : ( tn_-1 ) ) {{ q41_u_code_i[ i ] = 10000; q45_u_code_i[ i ] = 50;   }}
      
      for ( i in 1                 : n_t_wgh_1 ) {{ q51_u_code_i[ i ] = 5; q61_u_code_i[ i ] = 1; }}
      for ( i in ( n_t_wgh_1 + 1 ) : ( tn_-1 ) ) {{ q51_u_code_i[ i ] = 5; q61_u_code_i[ i ] = 1; }}
      
      //// q14 flows
                                                                         sum_q1j               = b_rate_rel_to_blood_naive_ * brf_free_param;
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q14_min_I_i[ i ]      = ( 0                                                     ) / al_[1]; }}
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q14_max_I_i[ i ]      = ( al_[of+1] * ( q41_u_code_i[ i ] + q45_u_code_i[ i ] ) ) / al_[1]; }}
                                                                         q14_sum_min           = ( al_[1] * sum_q1j - al_[1] * q12_u_code ) / al_[1] ;
                                                                         q14_sum_max           = ( al_[1] * sum_q1j - al_[1] * 0          ) / al_[1] ;
      ////  FIRST FLOW
                                                                         q14_propos_i[ 1 ]     = q14_min_I_i[ 1 ] + ( min( { q14_max_I_i[ 1 ], q14_sum_max } ) - sum( q14_min_I_i )                  - 0                           ) * q14_base_i[ 1 ];
                                                                         q14_i[ 1 ]            = max( { 0, min( { q14_propos_i[ 1 ], q14_max_I_i[ 1 ] } ) } );
      if ( ( tn_-2 ) >= 2 ) {{ for ( i in 2 : ( tn_-2 ) ) {{ of = i * 3; q14_propos_i[ i ]     = q14_min_I_i[ i ] + ( min( { q14_max_I_i[ i ], q14_sum_max } ) - sum( q14_min_I_i[ i : ( tn_-1 ) ] ) - sum( q14_i[ 1 : ( i-1 ) ] ) ) * q14_base_i[ i ]; }} }}
      if ( ( tn_-2 ) >= 2 ) {{ for ( i in 2 : ( tn_-2 ) ) {{ of = i * 3; q14_i[ i ]            = max( { 0, min( { q14_propos_i[ i ], q14_max_I_i[ i ] } ) } ); }}  }}
                                                                         q14_min_last          = q14_sum_min - sum( q14_i[ 1 : ( tn_-2 ) ] );
      if ( q14_min_last >= 0 ) {{                                        q14_propos_i[ tn_-1 ] = q14_min_last + ( q14_sum_max - q14_sum_min ) * q14_base_i[ tn_-1 ];    }}
      if ( q14_min_last <  0 ) {{                                        q14_propos_i[ tn_-1 ] = ( q14_sum_max - sum( q14_i[ 1 : ( tn_-2 ) ] ) ) * q14_base_i[ tn_-1 ]; }}
                                                                         q14_i[ tn_-1 ]        = max( { 0, min( { q14_propos_i[ tn_-1 ], q14_max_I_i[ tn_-1 ] } ) } );
      //// END of q14 flows
                                                                         q12_tmp               = sum_q1j - sum( q14_i );
     
      
      //// q36 flows
                                                                         sum_al_6_q63          = 0;
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          sum_al_6_q63          = sum_al_6_q63 + al_[of+3] * q63_i[ i ]; }}
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q66_min_i[ i ]        = - q61_u_code_i[ i ] - q63_i[ i ] - q56_q65_u_code; }}
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q66_max_i[ i ]        = - 0                 - q63_i[ i ] - 0             ; }}
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q56_i[ i ]            = q56_base_i[ i ] * ( ( al_[of+3] / al_[of+2] ) * q63_i[ i ] ); }}
      
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q36_max_I_i[ i ]      =           ( al_[of+3] * ( -q66_min_i[ i ] ) - ( al_[of+2] * q56_i[ i ] ) ) / al_[3]; }}
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q36_min_I_i[ i ]      = max( { 0, ( al_[of+3] * ( -q66_max_i[ i ] ) - ( al_[of+2] * q56_i[ i ] ) ) / al_[3] } ); }}  
                                                                         q36_sum_min           = ( sum_al_6_q63 - al_[3] * q32 + al_[2] * 0                                      ) / al_[3] ;
                                                                         q36_sum_max           = q36_sum_min                 + ( al_[2] * min( { q23_u_code, max_flow_to_N3_ } ) ) / al_[3] ;
      ////  FIRST FLOW
                                                                         q36_propos_i[ 1 ]     = q36_min_I_i[ 1 ] + max( { 0, ( min( { q36_sum_max, q36_max_I_i[ 1 ] } ) - sum( q36_min_I_i ) ) } ) * q36_base_i[ 1 ];
                                                                         q36_i[ 1 ]            = max( { 0, min( { q36_propos_i[ 1 ], q36_max_I_i[ 1 ] } ) } );
      ////  SECOND TO LAST-BUT-ONE FLOW - if present
      if ( ( tn_-2 ) >= 2 ) {{ for ( i in 2 : ( tn_-2 ) ) {{ of = i * 3; q36_propos_i[ i ]     = q36_min_I_i[ i ] + ( min( { q36_sum_max, q36_max_I_i[ i ] } ) - sum( q36_min_I_i ) - sum( q36_i[ 1 : ( i-1 ) ] ) ) * q36_base_i[ i ]; }} }}
      if ( ( tn_-2 ) >= 2 ) {{ for ( i in 2 : ( tn_-2 ) ) {{ of = i * 3; q36_i[ i ]            = max( { 0, min( { q36_propos_i[ i ], q36_max_I_i[ i ] } ) } ); }}  }}
      ////  LAST FLOW
                                                                         q36_min_last          = q36_sum_min - sum( q36_i[ 1 : ( tn_-2 ) ] );
                                                                         q36_max_last          = q36_sum_max - sum( q36_i[ 1 : ( tn_-2 ) ] );
      if ( q36_min_last >= 0 ) {{                                        q36_propos_i[ tn_-1 ] = q36_min_last + ( min( { q36_max_last, q36_max_I_i[ tn_-1 ] } ) - q36_min_last ) * q36_base_i[ tn_-1 ]; }}
      if ( q36_min_last <  0 ) {{                                        q36_propos_i[ tn_-1 ] = q36_max_last * q36_base_i[ tn_-1 ]; }}
                                                                         q36_i[ tn_-1 ]        = max( { 0, min( { q36_propos_i[ tn_-1 ], q36_max_I_i[ tn_-1 ] } ) } );
      //// END of q36 flows
      

      //// flows going to and from N3 and N6:
                                                                         q23                   = max( { 0, ( al_[3] * ( q32 + sum( q36_i ) ) - sum_al_6_q63 ) / al_[2] } );
      if ( max_flow_to_N3_ < 1 ) {{                                      q39                   = q39 * max_flow_to_N3_; }}
    
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                      q61_q65_sum_propos_i[ i ] = ( al_[3] * q36_i[ i ] + al_[of+2] * q56_i[ i ] - al_[of+3] * q63_i[ i ]) / al_[of+3]; }}    
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q61_q65_sum_i[ i ]    = max( { 0, q61_q65_sum_propos_i[ i ] } ); }}
      
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q65_i[ i ]            = min( { q56_q65_u_code, q65_base_i[ i ] * ( q61_q65_sum_i[ i ] ) } ); }}
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q61_i[ i ]            = q61_q65_sum_i[ i ] - q65_i[ i ]; }}
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q66_i[ i ]            = - q61_i[ i ] - q63_i[ i ] - q65_i[ i ]; }}
      //// END of flows going to and from N3 and N6:
      
      
      //// q25 flows
                                                                         q12_u_code = q12_tmp; q12_l_code = q12_tmp; 
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q45_max_i[ i ]        = max( { 0, ( al_[1] * q14_i[ i ] ) / al_[of+1] } ); }}
      
       
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q25_min_N5_i[ i ]     = max( { 0, al_[of+2] * ( 0                 + q56_i[ i ] + 0          ) - al_[of+3] * q65_i[ i ] - al_[of+1] * q45_max_i[ i ] } ) / al_[2]; }}
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q25_max_N5_i[ i ]     =    (      al_[of+2] * ( q51_u_code_i[ i ] + q56_i[ i ] + q52_u_code ) - al_[of+3] * q65_i[ i ] - al_[of+1] * 0                ) / al_[2]; }}
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q25_max_N5_i[ i ]     = min( { q25_max_N5_i[ i ], 0.1 / al_[2] }  ); }}
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q25_i[ i ]            = q25_min_N5_i[ i ] + q25_base_i[ i ] * ( q25_max_N5_i[ i ] - q25_min_N5_i[ i ] ); }}
      
      al2_out = al_[2] * ( q23 + sum( q25_i ) );
      
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q52_min_N5_i[ i ]     = max( { 0,                 ( al_[of+3] * q65_i[ i ] + al_[of+1] * 0              + al_[2] * q25_i[ i ] - al_[of+2] * ( q51_u_code_i[ i ] + q56_i[ i ] ) ) / al_[of+2] } ); }}
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q52_max_N5_i[ i ]     = max( { q52_min_N5_i[ i ], ( al_[of+3] * q65_i[ i ] + al_[of+1] * q45_max_i[ i ] + al_[2] * q25_i[ i ] - al_[of+2] * ( 0                 + q56_i[ i ] ) ) / al_[of+2] } ); }}
      
      
      
                                                                         sum_al_5_q52_max_N5   = 0;
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          sum_al_5_q52_max_N5   = sum_al_5_q52_max_N5 + al_[ 2+3*i ] * q52_max_N5_i[ i ]; }}
                                                                         sum_al_5_q52_min_N5   = 0;
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          sum_al_5_q52_min_N5   = sum_al_5_q52_min_N5 + al_[ 2+3*i ] * q52_min_N5_i[ i ]; }}
      
      for ( i in 1 : (     1 ) ) {{ of = i * 3;                          q52_min_N2_i[ i ]     = max( { 0,                 ( al_[ 2+3*i ] * q52_max_N5_i[ i ] + al2_out - al_[1] * q12_u_code - al_[3] * q32 - sum_al_5_q52_max_N5 ) / al_[of+2] } ); }}
      for ( i in 1 : (     1 ) ) {{ of = i * 3;                          q52_max_N2_i[ i ]     = max( { q52_min_N2_i[ i ], ( al_[ 2+3*i ] * q52_min_N5_i[ i ] + al2_out - al_[1] * q12_l_code - al_[3] * q32 - sum_al_5_q52_min_N5 ) / al_[of+2] } ); }}
      for ( i in 1 : (     1 ) ) {{ of = i * 3;                          q52_min_i[ i ]        = max( { q52_min_N2_i[ i ], q52_min_N5_i[ i ] } ); }}
      for ( i in 1 : (     1 ) ) {{ of = i * 3;                          q52_max_i[ i ]        = min( { q52_max_N2_i[ i ], q52_max_N5_i[ i ] } ); }}
      for ( i in 1 : (     1 ) ) {{ of = i * 3;                          q52_i[ i ]            = q52_min_i[ i ] + q52_base_i[ i ] * ( q52_max_i[ i ] - q52_min_i[ i ] );  q52_min_N5_i[ i ] = q52_i[ i ];  q52_max_N5_i[ i ] = q52_i[ i ];  }}
      
      
      ////  SECOND TO LAST-BUT-ONE FLOW - if present
                                                                         sum_al_5_q52_max_N5   = 0;
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          sum_al_5_q52_max_N5   = sum_al_5_q52_max_N5 + al_[ 2+3*i ] * q52_max_N5_i[ i ]; }}
                                                                         sum_al_5_q52_min_N5   = 0;
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          sum_al_5_q52_min_N5   = sum_al_5_q52_min_N5 + al_[ 2+3*i ] * q52_min_N5_i[ i ]; }}
      
      if ( ( tn_-2 ) >= 2 ) {{ for ( i in 2 : ( tn_-2 ) ) {{ of = i * 3; q52_min_N2_i[ i ]     = max( { 0,                 ( al_[ 2+3*i ] * q52_max_N5_i[ i ] + al2_out - al_[1] * q12_u_code - al_[3] * q32 - sum_al_5_q52_max_N5 ) / al_[of+2] } ); }} }}
      if ( ( tn_-2 ) >= 2 ) {{ for ( i in 2 : ( tn_-2 ) ) {{ of = i * 3; q52_max_N2_i[ i ]     = max( { q52_min_N2_i[ i ], ( al_[ 2+3*i ] * q52_min_N5_i[ i ] + al2_out - al_[1] * q12_l_code - al_[3] * q32 - sum_al_5_q52_min_N5 ) / al_[of+2] } ); }} }}
      if ( ( tn_-2 ) >= 2 ) {{ for ( i in 2 : ( tn_-2 ) ) {{ of = i * 3; q52_min_i[ i ]        = max( { q52_min_N2_i[ i ], q52_min_N5_i[ i ] } ); }} }}
      if ( ( tn_-2 ) >= 2 ) {{ for ( i in 2 : ( tn_-2 ) ) {{ of = i * 3; q52_max_i[ i ]        = min( { q52_max_N2_i[ i ], q52_max_N5_i[ i ] } ); }} }}
      if ( ( tn_-2 ) >= 2 ) {{ for ( i in 2 : ( tn_-2 ) ) {{ of = i * 3; q52_i[ i ]            = q52_min_i[ i ] + q52_base_i[ i ] * ( q52_max_i[ i ] - q52_min_i[ i ] );  q52_min_N5_i[ i ] = q52_i[ i ];  q52_max_N5_i[ i ] = q52_i[ i ];  }} }}
      
      
      //// i = tn_-1
                                                                         sum_al_5_q52_max_N5   = 0;
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          sum_al_5_q52_max_N5   = sum_al_5_q52_max_N5 + al_[ 2+3*i ] * q52_max_N5_i[ i ]; }}
                                                                         sum_al_5_q52_min_N5   = 0;
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          sum_al_5_q52_min_N5   = sum_al_5_q52_min_N5 + al_[ 2+3*i ] * q52_min_N5_i[ i ]; }}
      
      for ( i in ( tn_-1 ) : ( tn_-1 ) ) {{ of = i * 3;                  q52_min_N2_i[ i ]     = max( { 0,                 ( al_[ 2+3*i ] * q52_max_N5_i[ i ] + al2_out - al_[1] * q12_u_code - al_[3] * q32 - sum_al_5_q52_max_N5 ) / al_[of+2] } ); }}
      for ( i in ( tn_-1 ) : ( tn_-1 ) ) {{ of = i * 3;                  q52_max_N2_i[ i ]     = max( { q52_min_N2_i[ i ], ( al_[ 2+3*i ] * q52_min_N5_i[ i ] + al2_out - al_[1] * q12_l_code - al_[3] * q32 - sum_al_5_q52_min_N5 ) / al_[of+2] } ); }}
      for ( i in ( tn_-1 ) : ( tn_-1 ) ) {{ of = i * 3;                  q52_min_i[ i ]        = max( { q52_min_N2_i[ i ], q52_min_N5_i[ i ] } ); }}
      for ( i in ( tn_-1 ) : ( tn_-1 ) ) {{ of = i * 3;                  q52_max_i[ i ]        = min( { q52_max_N2_i[ i ], q52_max_N5_i[ i ] } ); }}
      for ( i in ( tn_-1 ) : ( tn_-1 ) ) {{ of = i * 3;                  q52_i[ i ]            = q52_min_i[ i ] + q52_base_i[ i ] * ( q52_max_i[ i ] - q52_min_i[ i ] );  q52_min_N5_i[ i ] = q52_i[ i ];  q52_max_N5_i[ i ] = q52_i[ i ];  }}
      //// END of q25 flows
      
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q51_min_i[ i ]        = min( { q51_u_code_i[ i ], max( { 0,              ( al_[2] * q25_i[ i ] + al_[of+1] * 0              + al_[of+3] * q65_i[ i ] - al_[of+2] * ( q52_i[ i ] + q56_i[ i ] ) ) / al_[of+2]} ) } ); }}  
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q51_max_i[ i ]        = min( { q51_u_code_i[ i ], max( { q51_min_i[ i ], ( al_[2] * q25_i[ i ] + al_[of+1] * q45_max_i[ i ] + al_[of+3] * q65_i[ i ] - al_[of+2] * ( q52_i[ i ] + q56_i[ i ] ) ) / al_[of+2] } ) } ); }} 
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q51_i[ i ]            = q51_min_i[ i ] + q51_base_i[ i ] * ( q51_max_i[ i ] - q51_min_i[ i ] ); }}
      //// The rest is determined by previously sampled parameters:
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q55_i[ i ]            = - q51_i[ i ] - q52_i[ i ] - q56_i[ i ]; }}

      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q45_i[ i ]            = ( al_[of+2] * ( q51_i[ i ] + q52_i[ i ] + q56_i[ i ] ) - ( al_[2] * q25_i[ i ] + al_[of+3] * q65_i[ i ] ) ) / al_[of+1]; }}
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q14_i[ i ]            = max( { q14_i[ i ], q45_i[ i ] * al_[of+1] / al_[1] } ); }}     
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q41_i[ i ]            = max( { 0, ( al_[1] * q14_i[ i ] / al_[of+1] ) - q45_i[ i ] } ); }}
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          q44_i[ i ]            = - q41_i[ i ] - q45_i[ i ]; }}
      
                                                                         sum_al_5_q52          = 0;
      for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3;                          sum_al_5_q52          = sum_al_5_q52 + ( al_[of+2] * q52_i[ i ] ); }}
      q12 = max( { 0,  ( al_[2] * q23 + al_[2] * sum( q25_i ) - al_[3] * q32 - sum_al_5_q52 ) / al_[1] } );
    }}//REMOVE IN R
    ////  hack for the R parser    
                                                                         q33                   = - q32 - sum( q36_i );
                                                                         q22                   = - q23 - sum( q25_i );
                                                                         q11                   = - q12 - sum( q14_i );
      
                                                                         q11                   = q11 - q17;
                                                                         q22                   = q22 - q28;
                                                                         q33                   = q33 - q39;
  //// week 1 nuisance parameters are not estimated, rather forced to high value
  //// Please kaeep no spaces in af[i,j]
  for ( i in 1 : wn_ ) {{ for ( j in 1 : tn_ ) {{ af[i,j] = af0[i,j]; }} }}
  af[ 1, 2 ] = 999;
  af[ 1, 3 ] = 999;
  
  //// !!! DO NOT REWRITE THE FOLLOWING LINE - IT IS DETECTED BY R PARSER
  // end of nonfree variables
}


model {
  // populations
  
  // 01    host blood naive
  // 02    host blood activ
  // 03    host blood cd69p
  
  // 04    host tissue_i naive
  // 05    host tissue_i activ
  // 06    host tissue_i cd69p 
  // etc.
  
  // ho+01    host blood naive                                                  // ho ... Host Offset
  // ho+02    host blood activ
  // ho+03    host blood cd69p
  
  // ho+04    host tissue_i naive
  // ho+05    host tissue_i activ
  // ho+06    host tissue_i cd69p
  // etc.
  
  matrix[ d, d ] Q = rep_matrix( 0, d, d );
  vector[ d ] a;
  array[ wn_, tn_ ] vector[ 2*pn_ ] alpha;
  int of = -1000;                                                               // of ... tissue OFfset
  int ho = 3 * tn_;                                                             // ho ... Host Offset
  
  // start of Q alignment - qij values into matrix Q
  // donor animal:
  // blood graph node flows and donor->host flows
  Q[ 1, { 1, 2,    ho+1 } ] = [ q11, q12,      q17 ];
  Q[ 2, {    2, 3, ho+2 } ] = [      q22, q23, q28 ];
  Q[ 3, {    2, 3, ho+3 } ] = [      q32, q33, q39 ];
  // blood flow to tissue_i
  for( i in 1 : ( tn_-1 ) ) {{ of = i * 3; Q[ 1, of+1           ] = q14_i[i] ; }}                  // my R parser: "];" -> ");" -> the space is crucial here: "q14_i[i] ;"
  for( i in 1 : ( tn_-1 ) ) {{ of = i * 3; Q[ 2,      of+2      ] = q25_i[i] ; }}
  for( i in 1 : ( tn_-1 ) ) {{ of = i * 3; Q[ 3,           of+3 ] = q36_i[i] ; }}
  // tissue i graph node flows
  
  for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3; Q[ of+1, { 1,       of+1, of+2      , ho+1 } ] = [ q41_i[i],                              q44_i[i], q45_i[i]          , 0 * q41_i[i] ]; }}
  for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3; Q[ of+2, { 1, 2,          of+2, of+3, ho+1 } ] = [ q51_i[i],     q52_i[i],                          q55_i[i], q56_i[i], 0 * q51_i[i] ]; }}
  for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3; Q[ of+3, { 1,    3,       of+2, of+3, ho+1 } ] = [ q61_i[i],               q63_i[i],                q65_i[i], q66_i[i], 0 * q61_i[i] ]; }}
  
  // host animal:
  // blood graph node flows and host->donor flows
  Q[ ho+1, { 1, ho+1, ho+2       } ] = [ q17, q11, q12      ];
  Q[ ho+2, { 2,       ho+2, ho+3 } ] = [ q28,      q22, q23 ];
  Q[ ho+3, { 3,       ho+2, ho+3 } ] = [ q39,      q32, q33 ];
  // blood flow to tissue_i
  for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3; Q[ ho+1, ho+of+1                   ] = q14_i[i] ; }}
  for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3; Q[ ho+2,           ho+of+2         ] = q25_i[i] ; }}    // my R parser: there needs to be a space behind: "q25_i[i] ;"
  for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3; Q[ ho+3,                   ho+of+3 ] = q36_i[i] ; }}
  // tissue i graph node flows
  for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3; Q[ ho+of+1, { 1,    ho+1,             ho+of+1, ho+of+2          } ] = [ 0 * q41_i[i],      ( q41_i[i] ),                       q44_i[i], q45_i[i]           ]; }}
  for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3; Q[ ho+of+2, { 1,    ho+1, ho+2,                ho+of+2, ho+of+3 } ] = [ 0 * q51_i[i],        q51_i[i],     q52_i[i],                     q55_i[i], q56_i[i] ]; }}
  for ( i in 1 : ( tn_-1 ) ) {{ of = i * 3; Q[ ho+of+3, { 1,    ho+1,       ho+3,          ho+of+2, ho+of+3 } ] = [ 0 * q61_i[i],        q61_i[i],               q63_i[i],           q65_i[i], q66_i[i] ]; }}
  // end of Q alignment
  

  for ( j in 1 : wn_ )
  {
    a = matrix_exp( 7 * ( w_[ j ] - w0i_ ) * Q )' * a0i_;                     
  
        for ( k in 1 : tn_ )
        {
            alpha[ j , k ] = append_row( a[ (k-1)*pn_ + 1 : k*pn_ ],
                a[ (tn_+k-1)*pn_ + 1 : (tn_+k)*pn_ ] );
            alpha[ j , k ] *= af[ j, k ] / sum( alpha[ j, k ] );
        }
    }

    
    for ( i in 1 : n_ )
    {
        if ( wi_[ i ] == 1 && ti_[ i ] == 1 ) {
            target += dirichlet_lpdf( x_[ i ] | alpha[ wi_[ i ], ti_[ i ] ] ) * 0;
        } else {
            target += dirichlet_lpdf( x_[ i ] | alpha[ wi_[ i ], ti_[ i ] ] ) * 1;
        }
    }
}

