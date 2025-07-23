# SET_PATHS.r

SET_PATHS <- function( PROJECT_PATH ) 
{
  CODE_PATH = sprintf( '%s/%s', PROJECT_PATH, 'script' );
  RAW_PATH = sprintf( '%s/%s', PROJECT_PATH, 'raw' );                                                 
  PROCESSED_PATH = sprintf( '%s/%s', PROJECT_PATH, 'processed' );             
  ANALYSIS_PATH = sprintf( '%s/%s', PROJECT_PATH, 'analyzed' );
  ANALYSIS_DONE_PATH = sprintf( '%s/%s', PROJECT_PATH, 'analyzed_ZENODO' );
  RESULTS_PATH = sprintf( '%s/%s', PROJECT_PATH, 'results' );
  PS <- list( PROJECT_PATH, CODE_PATH, RAW_PATH, 
              PROCESSED_PATH, ANALYSIS_PATH, RESULTS_PATH,
              ANALYSIS_DONE_PATH )
  names( PS ) <- c( 'PROJECT_PATH', 'CODE_PATH', 'RAW_PATH', 
                    'PROCESSED_PATH', 'ANALYSIS_PATH', 'RESULTS_PATH',
                    'ANALYSIS_DONE_PATH' )
  return( PS )
}
