args <- base::commandArgs(trailingOnly=TRUE)

if (base::length(args) != 1) {
  args_string <- base::paste(args, collapse=', ')
  stop(base::paste0('expected exactly 1 argument, but got: ', args_string))
}

packages <- base::list()
is_paired <- FALSE
is_darts <- FALSE
if (args[1] == 'paired') {
  is_paired <- TRUE
  packages <- base::list(
                        list(name='doParallel'),
                        list(name='foreach'),
                        list(name='iterators'),
                        list(name='nloptr')
                    )
} else if (args[1] == 'darts') {
  is_darts <- TRUE
  packages <- base::list(
                        list(name='doSNOW'),
                        list(name='getopt'),
                        list(name='ggplot2'),
                        list(name='mixtools'),
                        list(name='Rcpp')
                    )
} else {
  stop('expected either "paired" or "darts" as the only argument')
}

repos <- 'http://cran.us.r-project.org'

for (i in 1:length(packages)) {
  package <- packages[[i]]
  installed <- require(package$name, character.only=TRUE)
  if (installed) {
    next
  }

  install.packages(package$name, repos=repos)
  installed <- require(package$name, character.only=TRUE)
  if (!installed) {
    stop(paste('could not install', package$name))
  }
}

if (is_paired) {
  pairadise_installed <- require('PAIRADISE')
  if (!pairadise_installed) {
    install.packages('./PAIRADISE/pairadise/src/pairadise_model/',
                     repos=NULL)
  }
  pairadise_installed <- require('PAIRADISE')
  if (!pairadise_installed) {
    stop('could not install pairadise')
  }
} else if (is_darts) {
  darts_installed <- require('Darts')
  if (!darts_installed) {
    ## --preclean removes existing compilation results.
    ## This avoids an error if the downloaded DARTS code included build results.
    install.packages('./DARTS/Darts_BHT/Darts_BHT/Darts/',
                     repos=NULL,
                     INSTALL_opts=c('--preclean'))
  }
  darts_installed <- require('Darts')
  if (!darts_installed) {
    stop('could not install darts')
  }
}
