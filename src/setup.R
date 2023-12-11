repo <- "https://cloud.r-project.org"
# packrat::init()
# packrat::snapshot()
install.packages(c("tidyverse", "BiocManager", "RColorBrewer", "pheatmap"),
                 repos = repo)
BiocManager::install(c("DESeq2", "topGO"))
