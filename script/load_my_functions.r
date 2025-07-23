# load_my_functions.r

load_my_functions <- function( ps )
{
  my_fun_scipts <- c( 
    "get_parabiosis_const.r",
    
    "create_dirs.r",
    "get_path.r",
    "get_mcmc_pars.r",
    "get_stan_pars.r",
    "markov_premodel_prints_eps.r", 
    "prep_for_markov_model.r",
    "estimate_model_functions.r",
    "estimate_one_markov_model.r",
    "model_likelihood.r",
    "select_mcmc_chains.r",
    
    "get_HDI.r",
    "calculate_HDI.r",
    "write_sel_models_file.r",
    "flows_exit_comparison.r",
    
    "add_plot_to_figs_rda.r",
    "figure_flow_titles.r", 
    "preproc_fs_for_figures_2_4_3_5.r",
    "plot_dt_distributions.r",
    "sim_tissue_state_model.r",
    "plot_figures_2_4.r",
    "plot_figures_3_5.r",
    "plot_figure_combined_4.r",
    "group_flows.r",
    "write_dwell_time.r",
    "format_dwell_time_table.r",
    "visualize_results.r",
    "plot_16_comb_dwells.r",
    "plot_fig_7_comb_dwells.r",
    "plot_exit_types.r",
    "plot_influx.r",
    "plot_all_celltypes_combined_dwells.r",
    
    "converting_Q_and_flow_tables.r",
    "prep_data_figures_6.r",
    "merge_fig5_subfigures.r",
    "plot_all_tissue_trajectories.r",
    "plot_all_tissue_trajectories_helper.r" )   
  for ( fun_script in my_fun_scipts )
  {
    source( sprintf( "%s/%s", ps$CODE_PATH, fun_script ) )
  }  
}

ggexport_my <- function( ..., no_messages = TRUE )
{
  if ( no_messages ) {
    suppressMessages( {
      ggpubr::ggexport( ... )
    } )
  } else {
    ggpubr::ggexport( ... )
  } 
}
