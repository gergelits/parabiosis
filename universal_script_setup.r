# universal_script_setup.r

FILL_IN_MY_PROJECT_PATH <- ""
MY_PROJECT_PATH <- ifelse( nchar( FILL_IN_MY_PROJECT_PATH ) != 0,
                           FILL_IN_PROJECT_PATH, getwd() )
source( sprintf( "%s/SET_PATHS.r", MY_PROJECT_PATH ) )
ps <- PS <- SET_PATHS( PROJECT_PATH = MY_PROJECT_PATH )
source( sprintf( "%s/%s", ps$CODE_PATH, "load_packages.r" ) )
load_packages()
source( sprintf( "%s/%s", ps$CODE_PATH, "load_my_functions.r" ) )
load_my_functions( ps = ps )
