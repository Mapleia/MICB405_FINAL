library(stringr)
library(DESeq2)
library(readxl)
library(tidyverse)
library(dplyr)

cd_project_df <-read_excel("references/cd_project.xlsx") %>%
  dplyr::select(sra_id, col_name, type) %>%
  mutate(sample = col_name, SRR_ID = sra_id) %>%
  dplyr::select(SRR_ID, col_name, type)

get_all_counts <- function(dir_path) {
  files <- list.files(path=tab_file_dir)
  list_df <- vector(mode="list", length(files))
  for (i in 1:length(files)) {
    tab_file <- files[i]
    full_tab_file <- stringr::str_interp("${dir}/${file}",
                                         list(dir=tab_file_dir, file=tab_file)
    )
    srr_id <- str_split(tab_file, "_", simplify=TRUE)
    print(srr_id[3])
    cols <- c('gene_id', 'total', 'antisense', 'sense')
    df <- read_delim(full_tab_file,
                     col_names=cols,
                     col_types = c("c", "i", "i", "i", "c"),
                     delim="\t",
                     skip=4) %>%
      add_column(SRR_ID=srr_id[3])
    list_df[[i]] <- df
  }

  master_df <- do.call(rbind, list_df)

  master_typed_df <- left_join(master_df,
                               cd_project_df,
                               by="SRR_ID"
  ) %>%
    mutate(type = factor(type))

  return (master_typed_df)
}

tab_file_dir <- "inputs/Cd_exposure_STAR/to_zip/read_counts"
master_df <- get_all_counts(tab_file_dir) %>%
  dplyr::filter(type ==  "Ctrl" | type == "Cd_3days") %>%
  arrange(desc(type), col_name) %>%
  mutate(sense = as.integer(sense)) %>%
  dplyr::select(sense, gene_id, col_name, type)

shape_df_deseq2 <- function(df) {
  # make the column with the sample names
  # the column names for sense instead
  pivoted_df <- df %>%
    pivot_wider(names_from = col_name, values_from = sense)

  # make a matrix
  cd_frame <- pivoted_df %>% as.data.frame()
  rownames(cd_frame) <- cd_frame[,1]
  cd_frame<-cd_frame[,-c(1)]

  dat_matrix<- as.matrix(cd_frame)

  return (dat_matrix)
}

cd_matrix <- shape_df_deseq2(master_df)

cd_matrix_df <- tidyr::as_tibble(cd_matrix, rownames="gene_name") %>%
  mutate(across(
    Ctrl_1:Cd_3days_3,
    ~ as.integer(.x)
  )) %>%
  filter(
    if_any(
      Ctrl_1:Cd_3days_3,
      ~.x < 0
    )
  )

# create a metadata df
metadata <- master_df %>%
  select(col_name, type) %>%
  unique()

metadata <- data.frame(row.names = metadata$col_name,
             condition = metadata$type)


all(colnames(cd_matrix) == rownames(metadata))

dds_matrix <- DESeqDataSetFromMatrix(countData = cd_matrix, #matrix
                                     colData = metadata, #metadata file
                                     design = ~ condition)

dds_matrix$condition <- relevel(dds_matrix$condition, ref = "Ctrl")
levels(dds_matrix$condition)
dds <- DESeq(dds_matrix)
saveRDS(dds, "outputs/dds.rds")
