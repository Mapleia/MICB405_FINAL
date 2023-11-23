is_packrat_installed <- require("packrat")

if (!is_packrat_installed) {
  install.packages("packrat")
}

library("packrat")
packrat::restore()