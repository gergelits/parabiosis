# visualize_results.r

visualize_results <- function( m_tisNA = NULL )
{
  plot_figures_3_5( m_tisNA = m_tisNA )
  try( plot_figures_3_5_all( m_tisNA = m_tisNA, subgroup_short = "entry" ) )
  try( plot_figures_3_5_all( m_tisNA = m_tisNA, subgroup_short = "apoptosis" ) )
  
  # rare issue with annotation_logticks() at first run
  try( plot_figures_2_4( m_tisNA = m_tisNA ) ) 
  
  try( plot_figures_2_4( m_tisNA = m_tisNA, qij_I = 2L ) )
  n_tissues <- sub( ".*_t", "", m_tisNA$model_ver )
  if ( n_tissues == 2 ) {
    try( plot_figures_2_4( m_tisNA = m_tisNA, qij_I = 3L ) )
  }
  
  write_dwell_time( m_tisNA = m_tisNA )
  
  plot_dt_distributions( m_tisNA = m_tisNA, cellstate = "cd69p" )
  plot_dt_distributions( m_tisNA = m_tisNA, cellstate = "naive" )
  plot_dt_distributions( m_tisNA = m_tisNA, cellstate = "activ" )
  
  plot_aggr_tissuegroup_cellstate_barplot( m_tisNA = m_tisNA )
  
  # For 9-states flow figure
  write_aggregated_group_flows( m_tisNA = m_tisNA, write_csv = TRUE )
  
  get_aggr_tissuegroup_cellstate_areas( m_tisNA = m_tisNA, write_csv = TRUE )
  
  get_aggr_total_counts_areas( m_tisNA = m_tisNA, write_csv = TRUE )
  
  merge_fig5_subfigures( m_tisNA = m_tisNA )
  
  plot_all_tissue_trajectories( m_tisNA = m_tisNA )
}
