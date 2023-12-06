library(tidyverse)
library(topGO)
library(stringr)

cd_chlamy <- read_csv("outputs/chlamy_cd_results.csv") %>%
  mutate(up_down = as.factor(up_down))

top_cd_chlamy <- cd_chlamy  %>%
  filter(up_down == "UP")
bot_cd_chlamy <- cd_chlamy  %>%
  filter(up_down == "DOWN")

geneID2GO <- readMappings("outputs/chlamy_GOIDs.tsv")
geneUniverse <- names(geneID2GO)

# =============================================
down_gene_list <- factor(as.integer(geneUniverse %in% bot_cd_chlamy$gene_id))
names(down_gene_list) <- geneUniverse
down_GO_data <- new("topGOdata",
                    description = "Chalmy_low_cadmium",
                    ontology = "BP",
                    allGenes = down_gene_list,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)
down_result <- runTest(down_GO_data,
                       algorithm = "weight01",
                       statistic = "fisher")
# summarize the results
down_summary <- GenTable(down_GO_data,
                         weight01 = down_result,
                         orderBy = "down_result",
                         ranksOf = "down_result",
                         topNodes = 10)

down_GO_filtered <- down_summary %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  filter(weight01 <= 0.05) %>%
  head(n = 20) %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term))
# Now let's extract the order of the term column
down_order_term <- down_GO_filtered %>% 
  pull(Term) # pull() extracts a column as a vector
down_GO_filtered %>% 
  ggplot(aes(x= Term, y = GeneRatio, colour = weight01)) +
  geom_col(width = 0.05) +
  geom_point(size = 3) +
  coord_flip() +
  scale_x_discrete(limits = down_order_term) + 
  scale_colour_gradient(low = "red", high = "blue")

# =========================================

up_gene_list <- factor(as.integer(geneUniverse %in% top_cd_chlamy$gene_id))
names(up_gene_list) <- geneUniverse
up_GO_data <- new("topGOdata",
                  description = "Chalmy_high_cadmium",
                  ontology = "BP",
                  allGenes = up_gene_list,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO)
up_result <- runTest(up_GO_data,
                     algorithm = "weight01",
                     statistic = "fisher")

resultKS <- runTest(up_GO_data, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(up_GO_data, algorithm = "elim", statistic = "ks")

up_summary <- GenTable(up_GO_data,
                       weight01 = up_result,
                       orderBy = "up_result",
                       ranksOf = "up_result",
                       topNodes = 10)

up_GO_filtered <- up_summary %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  # filter(weight01 <= 0.05) %>%
  head(n = 20) %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term))

up_order_term <- up_GO_filtered %>% 
  pull(Term)

up_GO_filtered %>% 
  ggplot(aes(x= Term, y = GeneRatio, colour = weight01)) +
  geom_col(width = 0.05) +
  geom_point(size = 3) +
  coord_flip() +
  scale_x_discrete(limits = up_order_term) + 
  scale_colour_gradient(low = "red", high = "blue")
