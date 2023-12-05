library(readxl)
library(dplyr)
library(BiocManager)
library(reutils)

references_xlsx <- read_excel("references/references.xlsx", sheet = "list") %>%
  filter(`included?` == TRUE) 

for (project in references_xlsx$BioProject[1:1]) {
  sras <- esearch(project, db="sra") %>%
    efetch(db="sra", rettype="full")
  print(sras)
}

library(SRAdb)
sqlfile <-'SRAmetadb.sqlite'
if(!file.exists('SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile()
sra_con <- dbConnect(SQLite(),sqlfile)
