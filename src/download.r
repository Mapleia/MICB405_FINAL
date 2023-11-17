library(readxl)
library(dplyr)
library(BiocManager)

references_xlsx <- read_excel("references/references.xlsx", sheet = "list") %>%
  filter(`included?` == TRUE) 
