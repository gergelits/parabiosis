# create_dirs.r

create_dirs <- function( ps )
{
  for ( i in 1 : length( ps ) ) {
    if ( !dir.exists( ps[[ i ]] ) ) { 
      dir.create( ps[[ i ]], recursive = TRUE ) } 
  }
}
