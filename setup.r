is_packrat_installed <- require("packrat")

if (!is_packrat_installed) {
  install.packages("packrat")
}

library("packrat")
packrat::restore()

# This is required to run the notebooks using R-Kernel.
# PATH variable needs a a valid instance of jupyter for the below
# command to work.
# The setup of this is not covered in this repository as it
# differs from system to system.

# IRkernel::installspec()