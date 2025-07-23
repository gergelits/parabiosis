# write_sel_models_file.r

write_sel_models_file <- function( 
    m, celltype, ps, mcmc_pars, 
    model_name, model_ver, 
    sel_models_file_name )
{
  pmc <- get_pmc()
  dSelected_models <- tibble( 
    celltype = celltype,
    tissue = pmc$dTissueAllOrderedGroup$tissue.all.ordered, 
    model_id = NA, 
    model_dir_r = NA )
  
  for ( tissue in pmc$dTissueAllOrderedGroup$tissue.all.ordered )
  {
    m <- model( m = m, tissue = tissue )
    model_dir_r <- get_path( "model_dir_r", m = m, 
                             model_name = model_name, model_ver = model_ver,
                             mcmc_pars = mcmc_pars )
    model_dir <- get_path( 
      "model_dir", m = m, ps = ps, celltype = celltype, tissue = tissue,
      model_name = model_name, model_ver = model_ver,
      mcmc_pars = mcmc_pars )
    if (
      # this celltype-tissue-model is ready to be plotted
      file.exists(
        sprintf( "%s/parabio_fit%i.rda", model_dir, mcmc_pars$mcmc.iter.n ) ) &
      file.exists( 
        sprintf( "%s/parabio_fit_HDI.csv", model_dir ) )
    ) {
      dSelected_models$model_id[ dSelected_models$tissue == tissue ] <- 
        str_match( model_dir_r, "model_(.*?)_")[ , 2 ]
      dSelected_models$model_dir_r[ dSelected_models$tissue == tissue ] <- 
        model_dir_r
    }
  }
  dSelected_models %>% write.table( 
    sprintf( "%s/%s", ps$ANALYSIS_PATH, sel_models_file_name ), 
    sep = ",", row.names = FALSE )
}
 