library(DESeq2)
library(pheatmap)

dds <- readRDS("outputs/dds.rds")

vsd <- assay(vst(dds))
Z <- t(scale(t(vsd)))

topVarGenes <- rowVars(vsd) %>% 
  order(decreasing = TRUE) %>% 
  head(20) 

Z_topGenes <- Z[topVarGenes,] 
pheatmap(Z_topGenes,
         main = "Top 20 Genes")