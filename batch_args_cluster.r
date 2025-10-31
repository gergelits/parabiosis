# batch_args_cluster.r
if ( COMMENTS <- FALSE ) {
# Interface useful for model calling from shell script

# possible to set the paths of the project
# if not changed, the current path (getwd()) is used
}

FILL_IN_MY_PROJECT_PATH <- ""
MY_PROJECT_PATH <- ifelse( nchar( FILL_IN_MY_PROJECT_PATH ) != 0,
                         FILL_IN_PROJECT_PATH, getwd() )
source( sprintf( "%s/SET_PATHS.r", MY_PROJECT_PATH ) )
PS <- SET_PATHS( PROJECT_PATH = MY_PROJECT_PATH )

# set cell type and (typically) one tissue to estimate the Markov Model for
FF_CELLTYPE <- "Tconv"  
FF_TISSUE <- "Brain"    # tissue of interest
FF_TIS_0 <- "NULL"      # any extra compartment (good choice is also "Spleen")
FF_MODEL_ID <- "m0003"  # id of Markov chain model coded in Stan
FF_MCMC <- "i1001"       # file id with number of MCMC Markov chain simulations
FF_N_SIM <- 1001        # number of synthetic cells for simulating 
                        # cell behaviour based on already MCMC-estimated model

# optional parsing parameters from .sh script
args <- commandArgs( TRUE )
if ( length( args ) >= 1 ) { FF_CELLTYPE <- args[[ 1 ]] }
if ( length( args ) >= 2 ) { FF_TISSUE <- args[[ 2 ]] }
if ( length( args ) >= 3 ) { FF_TIS_0 <- args[[ 3 ]] }
if ( length( args ) >= 4 ) { FF_MCMC <- args[[ 4 ]] }
if ( length( args ) == 5 ) { FF_MODEL_ID <- args[[ 5 ]] }

# finding the .stan file based on its id
FF_MODEL_NAME <- sub( 
  pattern = "\\.stan", replacement = "", 
  list.files( path = PS$CODE_PATH, 
              pattern = sprintf( "model_%s.*\\.stan", FF_MODEL_ID ) ) )

# run the calculations:
source( sprintf( "%s/%s", PS$CODE_PATH, "calculate_parabiosis_paper.r" ) )
calculate_parabiosis_paper(
  ps = PS,
  CELLTYPE = FF_CELLTYPE,
  F_TISSUES = FF_TISSUE,
  VECTOR_TISSUES_0 = FF_TIS_0,
  MODEL_NAME = FF_MODEL_NAME,
  stan_pars_v = "v01",
  mcmc_pars_v = FF_MCMC,
  n_sim = FF_N_SIM )
