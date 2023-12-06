library(readr)
library(ggplot2)

res_filtered_final <- read_csv("outputs/chlamy_cd_results.csv")

res_filtered_final %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = up_down)) + 
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed")
