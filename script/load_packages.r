# load_packages.r

load_packages <- function()
{
  suppressPackageStartupMessages( {
    library( digest )
    library( expm )
    library( gridExtra )   
    library( gtools )      # gtools used also for ddirichlet()
    library( modeest )
    library( rstan )
    library( stringr )     
    library( scales )
    # library( tidyverse )     
    # Should call the individual packages instead.
    # https://www.tidyverse.org/blog/2018/06/tidyverse-not-for-packages/
    library( magrittr )
    library( ggplot2 )
    library( dplyr )
    library( readr )
    library( tidyr )
    library( forcats )
    library( patchwork )
    library( ggh4x )
  } )
  # can be used depending on the tidyverse / readr library version loaded
  try( options( readr.show_col_types = FALSE ) )  
  
  # the package posterior needed but better to avoid loading its functions 
  # and masking the default functions ( e.g., mad() )
  # similar for ggpubr
  stopifnot( nchar( system.file( package = "posterior" ) ) > 0 )
  stopifnot( nchar( system.file( package = "ggpubr" ) ) > 0 )
}
