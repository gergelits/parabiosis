# preproc_fs_for_figures_2_4_3_5.r

get_dAllTissues_Models_HDI_simple <- function( 
    ps, celltype, selected_models_file )
{
  if ( ! file.exists( selected_models_file ) ) { return( NULL ) }
  selected_models_table <- read_load_multi( fun = "read.csv", 
                                            full_path = selected_models_file )
  
  dAllTissues_Models_HDI <- NULL
  for ( i in 1 : nrow( selected_models_table ) ) {
    i.tissue <- selected_models_table$tissue[ i ]
    i.model_dir_r <- selected_models_table$model_dir_r[ i ]
    if ( !is.na( i.model_dir_r ) ) {
      dAllTissues_Models_HDI <- 
        dAllTissues_Models_HDI %>%
        bind_rows( read_csv( sprintf( "%s/%s/%s/%s/parabio_fit_HDI.csv",
                                      ps$ANALYSIS_PATH, celltype, 
                                      i.tissue, i.model_dir_r ) ) %>%
                     mutate( celltype = celltype,
                             tissue = i.tissue,
                             model_dir_r = i.model_dir_r ) %>%
                     dplyr::select( celltype, tissue, model_dir_r, 
                                    everything() ) ) 
    } else { 
      # print( sprintf( "WARNING: The tissue %s is missing", i.tissue ) ) 
    }
  }
  return( dAllTissues_Models_HDI )
}
