# group_flows.r

read_and_clean_total_counts_data <- function( 
    ps, celltype, count_type = "Median" )
{
  dTotal.counts <- 
    read_csv( sprintf( "%s/Total_counts/parabiosis_model_input_%s_counts.csv",
                     ps$PROCESSED_PATH, celltype ) ) %>% 
    dplyr::rename( tissue = "Tissue",
                   total.count = count_type,
                   total.count.sd = `Standard deviation` )
  return( dTotal.counts )
}


get_flows_all_tissues <- function( m_tisNA, ps, celltype, flow.part,
                                   model_name, model_ver, 
                                   mcmc_pars )
{
  if ( flow.part == "exit" ) { 
    get.flow.table <- get_exit_flow_table } else {
      get.flow.table <- get_entry_flow_table }
  tissues <- get_pmc()$TISSUES.NO.GUT
  model_dir_r <- get_path( "model_dir_r", m = m_tisNA,
                           model_name = model_name, model_ver = model_ver,
                           mcmc_pars = mcmc_pars )
  flow.table.all.tissues <- NULL
  for ( i.tissue in tissues )
  {  
    i.tissue.flow.table <- 
      get.flow.table( ps = ps, celltype = celltype, 
                      tissue = i.tissue, model_dir_r = model_dir_r ) 
    if ( !is.null( i.tissue.flow.table ) ) {
      i.tissue.flow.table <- i.tissue.flow.table %>% 
        mutate( tissue = i.tissue )
    }
    flow.table.all.tissues <- flow.table.all.tissues %>% 
      bind_rows( ., i.tissue.flow.table )
  }
  return( flow.table.all.tissues )
}


get_tissue_average_flows <- function( m_tisNA, ps, celltype, flow.part,
                                      model_name, model_ver, 
                                      mcmc_pars )
{
  flow.table <- get_flows_all_tissues( 
    m_tisNA = m_tisNA, ps = ps, celltype = celltype, flow.part = flow.part,
    model_name = model_name, model_ver = model_ver, 
    mcmc_pars = mcmc_pars )
  total.counts.table <- read_and_clean_total_counts_data(                      
    ps = ps, celltype = celltype )
  
  if ( all( dim( flow.table ) == c ( 0, 0 ) ) ) {
    return( tibble() )
  }
  
  tmp <- 
    flow.table %>% 
    left_join( ., get_pmc()$dTissueAllOrderedGroup.f %>%                                 
                 dplyr::select( tissue.all, tissue.group, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) ) %>%
    left_join( ., total.counts.table %>% dplyr::select( tissue, total.count ), 
               by = "tissue" )
  
  average.flows.table <-
    tmp %>% 
    # par used again technical         
    group_by( f.tissue.group, par ) %>%                                        
    mutate( sum.total.count.per.tissue.group = sum( total.count ) ) %>% 
    mutate( weight.of.tissue.in.tissue.group = 
              total.count / sum.total.count.per.tissue.group ) %>% 
    mutate( flow.mean.w.by.tissue.in.tissue.group = 
              sum( weight.of.tissue.in.tissue.group * flow ) ) %>%
    
    mutate( median.flow = median( flow ) ) %>% 
    mutate( mean.flow = mean( flow ) ) %>% 
    mutate( log10.median.flow = log10( median.flow ) ) %>% 
    mutate( log10.median.flow.1000 = log10( median.flow * 1000 ) ) %>% 
    
    arrange( par, f.tissue.group ) %>% 
    ungroup() %>% 
    dplyr::select( f.tissue.group, population.1, population.2, par, 
                   mean.flow, flow.mean.w.by.tissue.in.tissue.group, flow.part, 
                   median.flow, log10.median.flow, log10.median.flow.1000 ) %>% 
    distinct() %>% 
    
    group_by( population.1, f.tissue.group ) %>% 
    mutate( sum.flow.mean.w.per.population.1 = 
              sum( flow.mean.w.by.tissue.in.tissue.group ) ) %>% 
    ungroup() %>% 
    dplyr::select( f.tissue.group, population.1, population.2, par, 
                   mean.flow, flow.mean.w.by.tissue.in.tissue.group, flow.part, 
                   sum.flow.mean.w.per.population.1, 
                   median.flow, log10.median.flow, log10.median.flow.1000 )
    
  return( average.flows.table )
}
