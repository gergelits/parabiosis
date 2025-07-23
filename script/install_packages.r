# install_packages.r

packages <- c( 
  "digest", "expm", "modeest", "rstan", "posterior",
  "ggplot2", "dplyr", "readr", "tidyr", "magrittr", "stringr",
  "ggpubr", "gridExtra", "gtools", "scales", "patchwork", "ggh4x" )

for ( package_i in packages ) {
  if ( ! requireNamespace( package_i, quietly = TRUE ) ) {
    install.packages( package_i, repos = "https://cloud.r-project.org/" )
  }
}
