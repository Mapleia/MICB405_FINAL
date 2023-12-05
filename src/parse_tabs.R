library(stringr)
library(DESeq2)
library(readxl)
library(tidyverse)

cd_project_df <-read_excel("references/cd_project.xlsx") %>%
  dplyr::select(sra_id, col_name, type) %>%
  rename(col_name = "sample") %>%
  rename(sra_id = "SRR_ID")




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

tab_file_dir <- "./Cd_exposure_STAR/to_zip/read_counts"
master_df <- get_all_counts(tab_file_dir)

shape_df_deseq2 <- function(input_df, metadata) {
  # make the column with the sample names
  # the column names for sense instead
  pivoted_df <- master_df %>%
    filter(type == "Cd_3days" | type == "Ctrl") %>%
    dplyr::select(sense, gene_id, sample) %>%
    pivot_wider(names_from = sample, values_from = sense) %>%
    as.data.frame()

  # set row.names
  rownames(pivoted_df) <- pivoted_df[,1]
  pivoted_df[,1] <- NULL

  # make a matrix
  dat_matrix<- as.matrix(pivoted_df)
  # check that they match
  colnames(dat_matrix) == rownames(metadata)

  return (dat_matrix)
}

cd_matrix <- shape_df_deseq2(master_df, metadata)

# create a metadata df
metadata <- data.frame(row.names = colnames(cd_matrix),
                       condition = c("1_hour", "1_hour", "12_hour", "12_hour")
)


dds_matrix <- DESeqDataSetFromMatrix(countData = cd_matrix, #matrix
                                     colData = metadata, #metadata file
                                     design = ~condition)

dds_matrix$condition <- relevel(dds_matrix$condition, ref = "Ctrl")
levels(dds_matrix$condition)
dds <- DESeq(dds_matrix)
saveRDS(dds, "outputs/dds.rds")

