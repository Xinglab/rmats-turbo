repos <- "http://cran.us.r-project.org"

packages <- list(
  list(name="nloptr"),
  list(name="foreach"),
  list(name="doParallel")
)

for (i in 1:length(packages)) {
  package <- packages[[i]]
  installed <- require(package$name, character.only=TRUE)
  if (installed) {
    next
  }

  install.packages(package$name, repos=repos)
  installed <- require(package$name, character.only=TRUE)
  if (!installed) {
    stop(paste("could not install", package$name))
  }
}

pairadise_installed <- require("PAIRADISE")
if (!pairadise_installed) {
  install.packages("./PAIRADISE/pairadise/src/pairadise_model/",
                   repos=NULL)
}
pairadise_installed <- require("PAIRADISE")
if (!pairadise_installed) {
  stop("could not install pairadise")
}
