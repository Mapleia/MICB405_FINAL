library(DESeq2)
library(RColorBrewer)
library(pheatmap)


dds <- readRDS("outputs/dds.rds")
rld <- rlog(dds)
plotPCA(rld, intgroup = "condition")

# Calculate distances between samples in our log-transformed data
sample_dists <- dist(t(assay(rld)))

# Convert the output to a matrix
sample_dist_matrix <- as.matrix(sample_dists)

# Remove the column names of our matrix
colnames(sample_dist_matrix) <- NULL

# Set the colour palette for our heatmap
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Generate a heatmap using the pheatmap package
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colours)
resultsNames(dds)
results(dds, contrast = c("condition", "Cd_3days"))

res <- results(dds, name = "condition_Cd_3days_vs_Ctrl") %>% as.data.frame()

res_no_NA <- res %>%
  drop_na()
glimpse(res_no_NA)

res_filtered <- res_no_NA %>%
  filter(padj <= 0.05)

glimpse(res_filtered)

res_filtered_final <- res_filtered %>%
  filter(log2FoldChange <= -1 | log2FoldChange >= 1) %>% # the '|' stand for OR here!
  rownames_to_column("gene_id") %>%
  mutate(up_down = case_when(
    padj < 0.05 & log2FoldChange > 2 ~ "UP", # label genes in the top right quadrant of the volcano plot
    padj < 0.05 & log2FoldChange < -2 ~ "DOWN", # label genes in the top left quadrant of the volcano plot
    TRUE ~ "NONE"
  ))

head(res_filtered_final)
write_csv(res_filtered_final, "outputs/chlamy_cd_results.csv")
