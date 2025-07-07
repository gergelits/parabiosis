# add_plot_to_figs_rda.r

add_plot_to_figs_rda <- function( m = NULL,
                                  gg_add, gg_add_name, 
                                  gg_figs_rda_name )
{
  full_gg_figs_rda_name <- sprintf( 
    "%s/%s", get_path( "figs_m_rda", m = m ), gg_figs_rda_name )
  
  if ( file.exists( full_gg_figs_rda_name ) ) {
    # NOTE: read_load_multi() cannot be used
    load( full_gg_figs_rda_name )
  } else {                                             
    # create new list, gg_figs_all
    gg_figs_all <- list()
  }
  if ( ! ( gg_add_name %in% names( gg_figs_all ) ) ) {  
    # add gg_add to gg_figs_all
    gg_figs_all[[ length( gg_figs_all ) + 1 ]] <- gg_add
    names( gg_figs_all )[ length( gg_figs_all ) ] <- gg_add_name
  } else {                                             
    # change gg_add in gg_figs_all
    # cat( sprintf( "Updating %s in %s\n", gg_add_name, gg_figs_rda_name ) )
    gg_figs_all[[ which( names( gg_figs_all ) == gg_add_name ) ]] <- gg_add
  }
  save( gg_figs_all, file = full_gg_figs_rda_name )
}
