# write_dwell_time.r 

write_dwell_time <- function( m_tisNA = NULL )
{
  celltype <- m_tisNA$celltype
  ps <- m_tisNA$ps
  mcmc_pars <- m_tisNA$mcmc_pars
  
  sel_models_file_name <- get_sel_models_file_name( m_tisNA )
  model_id <- get_model_id( sel_models_file_name )
  figs_dir <- get_path( "figs_m", m = m_tisNA )
  fig_filename_start <- get_fig_fname_start( m = m_tisNA )
  dwell_csv_full <- sprintf( 
    "%s/%s_dwell_CrI.csv", figs_dir, fig_filename_start )
  
  dAllTissues_Models_HDI <- get_dAllTissues_Models_HDI_simple( 
    ps = ps, celltype = celltype, selected_models_file = 
      sprintf( "%s/%s", ps$ANALYSIS_PATH, sel_models_file_name ) ) 
  
  dAllTissues_Models_HDI %>% 
    left_join( get_pmc()$dTissueAllOrderedGroup.f %>%                                  
                 dplyr::select( tissue.all, f.tissue, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) ) %>% 
    # multi-tissue models
    mutate( par_orig = par ) %>%                                                 
    mutate( par = gsub( pattern = "(q[0-9]{2,2})(_i\\[1\\])", replacement = "\\1", 
                        x = par_orig ) ) %>% 
    filter( par %in% c( "q44", "q55", "q66" ) ) %>% 
    mutate( dwell_time_mode_DAYS = -1 / mode,
            dwell_time_hdi_80_upper = -1 / hdi.80.lower,
            dwell_time_hdi_80_lower = -1 / hdi.80.upper,
            dwell_time_hdi_50_upper = -1 / hdi.50.lower,
            dwell_time_hdi_50_lower = -1 / hdi.50.upper ) %>% 
    mutate( cellstate = ifelse( 
      par == "q66", "CD69+",
      ifelse( par == "q55", "Antigen-experienced",
              ifelse( par == "q44", "Naive", "OTHER" ) ) ) ) %>% 
    dplyr::select( celltype, tissue, cellstate, 
                   starts_with( "dwell_time" ), par ) %>% 
    arrange( desc( par ), tissue ) %>% 
    dplyr::distinct() %>% 
    write_csv( ., dwell_csv_full ) 
  
    format_dwell_time_table( csv_name = dwell_csv_full, ps = ps )
    
    plot_dwell_times( m_tisNA = m_tisNA, 
                      csv_name = dwell_csv_full,
                      sel_models_file_name = sel_models_file_name )
}
