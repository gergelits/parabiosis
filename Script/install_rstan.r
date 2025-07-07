remove.packages( "rstan" )
library( rstan )


install.packages( "rstan" )

# downgrade?
if ( FALSE ) {
  # https://stackoverflow.com/questions/17082341/installing-older-version-of-r-package
packageurl <-
  "http://cran.r-project.org/src/contrib/Archive/rstan/rstan_2.26.23.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
}
# does not work.


# avoid V8?
if ( FALSE ) {
# https://github.com/stan-dev/rstan/issues/831
devtools::install_github("makoshark/rstan", ref="develop", subdir="rstan/rstan")
# compilation terminated.
# make: *** [/usr/lib/R/etc/Makeconf:204: chains.o] Error 1
# ERROR: compilation failed for package ‘rstan’
# * removing ‘/home/vaclavgergelits/R/x86_64-pc-linux-gnu-library/4.4/rstan’
# Warning message:
#   In i.p(...) :
#   installation of package ‘/tmp/RtmpNQnz8z/file2551587b3535/rstan_2.21.2.tar.gz’ had non-zero exit status
}
