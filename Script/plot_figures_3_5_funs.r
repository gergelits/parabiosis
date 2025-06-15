# plot_figures_3_5_funs.r

get.tests.table <- function( dd ) {
  dd <- 
    dd %>% mutate(
    id.tg.fs = factor( sprintf( "%s__%s", f.tissue.group, f.flow.subgroup ) ) ) %>% 
    filter( ! ( par %in% c( "q12", "q23", "q32" ) ) )
    
  
  d.wilcox.tests <- tidyr::crossing( id.tg.fs = levels( dd$id.tg.fs ) ) %>% 
    tidyr::separate( col = id.tg.fs, sep = "__", remove = FALSE, 
                     into = c( "tissue.group", "f.flow.subgroup" ) ) %>% 
    mutate( f.flow.subgroup = factor( f.flow.subgroup ) ) %>% 
    mutate( w.p.value = NA, w.stat = NA, 
            t.p.value = NA, t.stat = NA )
  
  for ( k in 1 : nlevels( dd$id.tg.fs ) ) {
    d.tg.fs <- 
      dd %>% 
      filter( id.tg.fs == levels( id.tg.fs )[ k ] ) %>% 
      mutate( f.cellstate.group = droplevels( f.cellstate.group ) )
      
    
    if ( nlevels( d.tg.fs$f.cellstate.group ) == 2 ) {
      w.t <- wilcox.test( 
        x = d.tg.fs$mode[ d.tg.fs$f.cellstate.group == 
                            levels( d.tg.fs$f.cellstate.group )[ 1 ] ],
        y = d.tg.fs$mode[ d.tg.fs$f.cellstate.group == 
                            levels( d.tg.fs$f.cellstate.group )[ 2 ] ], 
        paired = TRUE )
      k <- which( 
        d.wilcox.tests$id.tg.fs == levels( dd$id.tg.fs )[ k ] )
      
      d.wilcox.tests$w.p.value[ k ] <- w.t$p.value
      d.wilcox.tests$w.stat[ k ] <- w.t$statistic
    }
  }
  dd <- 
    dd %>% left_join( ., d.wilcox.tests, by = c( "id.tg.fs" = "id.tg.fs" ) )
    
  res <- list( dd, d.wilcox.tests )
  names( res ) <- c( "dd", "d.wilcox.tests" )
  return( res )
}


get.pval.annot.table <- function( d.stat.test, y.pos = 8e2 ) {
  if ( nrow( d.stat.test ) >= 2 ) {
    pval.annot.table <- 
      d.stat.test %>% 
      left_join( ., get_pmc()$dTissueAllOrderedGroup.f %>% 
                   dplyr::select( tissue.group, f.tissue.group.longer ), 
                 by = c( "tissue.group" = "tissue.group" ) ) %>% distinct() %>% 
      mutate( mode = y.pos ) %>% 
      filter( w.p.value <= 0.05 ) %>% 
      mutate( w.p.value.text = sprintf( "P=%s", signif( w.p.value, 2 ) ) )
      
    return( pval.annot.table ) } else
      return( NULL )
}
