# figure_flow_titles.r 

parabiosis_flow_titles <- function( celltype ) { 

  naive_resting <- ifelse( celltype == "Treg", "Resting", "Naive" )
  activ_antig_exp <- ifelse( celltype == "Treg", 
                             "Antigen-experienced", "Activated" )
  dFlowsTitles <-
  tibble(
    flow_between = c( "1_to_2", "1_to_4",
                      "2_to_3", "2_to_5", 
                      "3_to_2", "3_to_6",
                      "4_to_1", "4_to_5",
                      "5_to_1", "5_to_2", "5_to_6",
                      "6_to_1", "6_to_3", "6_to_5" ),
    from_or_to.titlestable = "to",
    flow.title = c( 
      sprintf( "%s activation in blood", celltype ), 
      sprintf( "%s %s entry", naive_resting, celltype ),
      sprintf( "blood.%s.activ -> blood.%s.cd69p (2 -> 3)", celltype, celltype ), 
      sprintf( "%s %s tissue entry", activ_antig_exp, celltype ),                          
      sprintf( "blood.%s.cd69p -> blood.%s.activ (3 -> 2)", celltype, celltype ), 
      sprintf( "CD69+ %s tissue entry", celltype ),
      sprintf( "%s %s leaving tissue or dying", naive_resting, celltype ),
      sprintf( "%s activation in tissue", celltype ),
      sprintf( "%s %s death rate", activ_antig_exp, celltype ), 
      sprintf( "%s %s leaving tissue", activ_antig_exp, celltype ), 
      sprintf( "%s gaining tissue residency phenotype (CD69+)", celltype ), 
      sprintf( "CD69+ %s death rate", celltype ), 
      sprintf( "CD69+ %s leaving tissue", celltype ),
      sprintf( "%s losing tissue residency (CD69+)", celltype ) ) ) %>% 
    
    bind_rows( ., tibble(
      flow_between = c( "1_from_4", "1_from_5", "1_from_6",
                        "2_from_1", "2_from_5", 
                        "4_from_1",
                        "5_from_2", "5_from_4", "5_from_6",
                        "6_from_5" ),
      from_or_to.titlestable = "from",
      flow.title = c( sprintf( "Arising %s naive (1 <- 4)", celltype ), 
                      sprintf( "Arising %s naive (1 <- 5)", celltype ), 
                      sprintf( "Arising %s naive (1 <- 6)", celltype ),
                      sprintf( "blood.%s.activ <- blood.%s.naive (2 <- 1)", celltype, celltype ), 
                      sprintf( "blood.%s.activ <- tissue.%s.activ (2 <- 5)", celltype, celltype ),
                      sprintf( "Naive %s entry (4 <- 1)", celltype ), 
                      sprintf( "Activated %s entry (5 <- 2)", celltype ), 
                      sprintf( "Naive %s activation (5 <- 4)", celltype ), 
                      sprintf( "tissue.%s.activ <- tissue.%s.CD69+ (5 <- 6)", celltype, celltype ), 
                      "Gain of residency phenotype (6 <- 5)" ) ) )
    
  
  dPopulationsTitles <- 
    tibble(
      population = c( 1 : 6 ),
      population.title = c( sprintf( "Blood naive %ss", celltype ), 
                            sprintf( "Blood activated %ss", celltype ), 
                            sprintf( "Blood CD69+ %ss", celltype ), 
                            sprintf( "Tissue naive %ss", celltype ), 
                            sprintf( "Tissue activated %ss", celltype ), 
                            sprintf( "Tissue CD69+ %ss", celltype ) ) )
    
  
  s.en <- "Tissue entry\n(rate rel. to blood counts)"
  s.ex <- "Tissue exit"
  s.a <- "Apoptosis"
  s.d <- "Differentation"
  s.diag <- "Diag"
  
  ss.dif <- "Differentation"
  ss.dedif <- "De-differentation"
  ss.a <- "Apoptosis"
  ss.ex <- "Tissue exit"
  ss.en <- "Tissue entry\n(rate rel. to blood counts)"
  ss.diag <- "Diag"
  
  c.n <- "Naive"
  c.a <- "Antigen-experienced"   
  c.r <- "Resident"
  c.diag <- "Diag"
  
  dFlowsGroups <-
    tibble( flow_between = c( "1_to_2", "1_to_4",
                              "2_to_3", "2_to_5", 
                              "3_to_2", "3_to_6",
                              "4_to_1", "4_to_5",
                              "5_to_1", "5_to_2", "5_to_6",
                              "6_to_1", "6_to_3", "6_to_5",
                              "q17", "q28", "q39", "brf_free_param", "sum_q1j" ), 
            par = c( "q12", "q14",
                     "q23", "q25",
                     "q32", "q36", 
                     "q41", "q45",
                     "q51", "q52", "q56",
                     "q61", "q63", "q65",
                     "q17", "q28", "q39", "brf_free_param", "sum_q1j" ),
            flow.group = c( s.d, s.en,
                            s.d, s.en,
                            s.d, s.en, 
                            s.ex, s.d,
                            s.a, s.ex, s.d,
                            s.a, s.ex, s.d,
                            s.diag, s.diag, s.diag, s.diag, s.diag ),
            f.flow.group = factor( flow.group, 
                                   levels = c( s.en, s.ex, s.a, s.d, s.diag ) ),
            flow.subgroup = c( ss.dif, ss.en,
                               ss.dif, ss.en,
                               ss.dedif, ss.en, 
                               ss.ex, ss.dif,
                               ss.a, ss.ex, ss.dif,
                               ss.a, ss.ex, ss.dedif,
                               ss.diag, ss.diag, ss.diag, ss.diag, ss.diag ),
            f.flow.subgroup = factor( 
              flow.subgroup, 
              levels = c( ss.en, ss.ex, ss.a, ss.dif, ss.dedif, ss.diag ) ),
            cellstate.group = c( c.n, c.n,   
                                 c.a, c.a,   
                                 c.r, c.r,   
                                 c.n, c.n,
                                 c.a, c.a, c.a,
                                 c.r, c.r, c.r,
                                 c.diag, c.diag, c.diag, c.diag, c.diag ),
            f.cellstate.group = factor( cellstate.group,
                                        levels = c( c.n, c.a, c.r, c.diag ) ) )
    
  
  dPopulationCodes <- tibble( 
    population.codes = c( 1:3 ),
    population = c( "host.treg.naive", "host.treg.activ", "host.treg.cd69p" ) )
    
  
  pft <- list( dFlowsTitles, dPopulationsTitles, dFlowsGroups, dPopulationCodes )
  names( pft ) <- c( "dFlowsTitles", "dPopulationsTitles", 
                     "dFlowsGroups", "dPopulationCodes" )
  return( pft )
}

get_gr <- function( type = NULL )
{
  if ( type == "cm_to_in" ) { return( 1 / 2.54 ) }
  
  if ( type == "f" ) { return( get_gr( "cm_to_in" ) * 1.7 ) }
 
  if ( type == "MIN_0" ) { return( 0.00001 ) }
  
  if ( type == "theme_blank" ) {
    f <- get_gr( "f" )
    theme_blank <- 
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
    return( theme_blank )
  }
  return( NULL )
}

# scales::hue()( 3 )[ 1 ]
cellstate_color <- c( "Naive" = scales::hue_pal()( 3 )[ 1 ],
                      "naive" = scales::hue_pal()( 3 )[ 1 ],
                      "resting" = scales::hue_pal()( 3 )[ 1 ], 
                      "Naive/Resting" = scales::hue_pal()( 3 )[ 1 ], 
                      "Antigen-experienced" = scales::hue_pal()( 3 )[ 2 ], 
                      "Resident" = scales::hue_pal()( 3 )[ 3 ],
                      "CD69+" = scales::hue_pal()( 3 )[ 3 ] )

cellstate_shape <- c( "Naive" = 19,               # circles
                      "naive" = 19,
                      "resting" = 19,
                      "Naive/Resting" = 19,
                      "Antigen-experienced" = 17, # triangles
                      "Resident" = 15,
                      "CD69+" = 15 )           # squares
