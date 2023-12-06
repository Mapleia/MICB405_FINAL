library(tidyverse)
library(topGO)
library(stringr)

# ===========================================
geneID2GO <- readMappings("outputs/chlamy_GOIDs.tsv")
geneUniverse <- names(geneID2GO)
# top, bottom order
topgo_parse <- function(top_df, bottom_df, names, limit) {

  up_gene_list <- factor(as.integer(geneUniverse %in% top_df$gene_id))
  down_gene_list <- factor(as.integer(geneUniverse %in% bottom_df$gene_id))
  
  names(up_gene_list) <- geneUniverse
  names(down_gene_list) <- geneUniverse
  
  # build the GOdata object in topGO for upregulated
  up_GO_data <- new("topGOdata",
                    description = names[1],
                    ontology = "BP",
                    allGenes = up_gene_list,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)
  
  # build the GOdata object in topGO for downregulated
  
  down_GO_data <- new("topGOdata",
                      description = names[2],
                      ontology = "BP",
                      allGenes = down_gene_list,
                      annot = annFUN.gene2GO,
                      gene2GO = geneID2GO)
  
  up_result <- runTest(up_GO_data,
                       algorithm = "weight01",
                       statistic = "fisher")
  down_result <- runTest(down_GO_data,
                         algorithm = "weight01",
                         statistic = "fisher")
  # summarize the results
  up_summary <- GenTable(up_GO_data,
                         weight01 = up_result,
                         orderBy = "up_result",
                         ranksOf = "up_result",
                         topNodes = limit)
  down_summary <- GenTable(down_GO_data,
                           weight01 = down_result,
                           orderBy = "down_result",
                           ranksOf = "down_result",
                           topNodes = limit)
  
  up_GO <- up_summary %>% 
    mutate(up_down = "UP")
  
  down_GO <- down_summary %>% 
    mutate(up_down = "DOWN")
  
  
  # Make a joined dataframe
  joined_GO_filtered_arranged <- bind_rows(up_GO, down_GO) %>%
    filter(weight01 <= 0.05) %>%
    mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
    arrange(GeneRatio) %>%
    mutate(Term = factor(Term)) %>%
    head(n = limit*2)
  
  return(joined_GO_filtered_arranged)
}

topgo_plot <- function(df) {
  # Extract the column order
  order_term_joined <- df %>% 
    pull(Term)
  
  joined_plot <- df %>% 
    ggplot(aes(x= Term, y = GeneRatio, color = weight01)) +
    geom_point(aes(size= Significant)) +
    coord_flip() +
    scale_x_discrete(limits = order_term_joined) +
    scale_color_gradient(low = "red", high = "blue") +
    theme_light() +
    labs(x = "GO Term Description", y = "Enrichment Ratio", color = "P-value", size = "Number of Significant Genes") +
    theme(panel.border = element_rect(color = "black"), 
          panel.grid = element_line(colour = "grey96"), 
          strip.background = element_rect(colour = "black")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
    facet_grid(.~ up_down)
  
  return(joined_plot)
}

# ===========================================
names_lead <- c("Chalmy_high_lead", "Chalmy_low_lead")
lead_chlamy <- read_csv("outputs/topGO/highlead_vs_ctrl.csv") %>%
  drop_na() %>%
  filter(padj <= 0.05) %>%
  filter(log2FoldChange <= -1 | log2FoldChange >= 1) %>% # the '|' stand for OR here!
  mutate(up_down = case_when(
    padj < 0.05 & log2FoldChange > 2 ~ "UP", # label genes in the top right quadrant of the volcano plot
    padj < 0.05 & log2FoldChange < -2 ~ "DOWN", # label genes in the top left quadrant of the volcano plot
    TRUE ~ "NONE"
  ))

top_lead_chlamy <- lead_chlamy  %>%
  filter(up_down == "UP")
bot_lead_chlamy <- lead_chlamy  %>%
  filter(up_down == "DOWN")

# bot_lead_chlamy <- read_csv("outputs/topGO/lowlead_vs_ctrl.csv")
lead_topgo <- topgo_parse(top_lead_chlamy, bot_lead_chlamy, names_lead, 10)
lead_plots <- topgo_plot(lead_topgo)
lead_plots

ggsave(
  "outputs/lead_topgo_joined_plot.png",
  lead_plots,
  height = 8,
  scale = 3,
  width = 10,
  units = "cm",
  dpi = 400,
)
# ===========================================
cd_chlamy <- read_csv("outputs/chlamy_cd_results.csv") %>%
  mutate(up_down = as.factor(up_down))

top_cd_chlamy <- cd_chlamy  %>%
  filter(up_down == "UP")
bot_cd_chlamy <- cd_chlamy  %>%
  filter(up_down == "DOWN")

names_cd <- c("Chalmy_Cd_up", "Chalmy_Cd_down")
cd_topgo <- topgo_parse(top_cd_chlamy, bot_cd_chlamy, names_cd, 10)
cd_plots <- topgo_plot(cd_topgo)

ggsave(
  "outputs/cad_topgo_joined_plot.png",
  cad_plots,
  height = 8,
  scale = 3,
  width = 10,
  units = "cm",
  dpi = 400,
)

# ===========================================

separate_plot <- function(down_summary, up_summary) {
  down_GO_filtered <- down_summary %>%
    mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
    filter(weight01 <= 0.05) %>%
    head(n = 20) %>% 
    arrange(GeneRatio) %>%
    mutate(Term = factor(Term))
  
  up_GO_filtered <- up_summary %>%
    mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
    filter(weight01 <= 0.05) %>%
    head(n = 20) %>% 
    arrange(GeneRatio) %>%
    mutate(Term = factor(Term))
  
  # Now let's extract the order of the term column
  down_order_term <- down_GO_filtered %>% 
    pull(Term) # pull() extracts a column as a vector
  
  up_order_term <- up_GO_filtered %>% 
    pull(Term)
  
  down_col_point_plot <- down_GO_filtered %>% 
    ggplot(aes(x= Term, y = GeneRatio, colour = weight01)) +
    geom_col(width = 0.05) +
    geom_point(size = 3) +
    coord_flip() +
    scale_x_discrete(limits = down_order_term) + 
    scale_colour_gradient(low = "red", high = "blue")
  
  up_col_point_plot <- up_GO_filtered %>% 
    ggplot(aes(x= Term, y = GeneRatio, colour = weight01)) +
    geom_col(width = 0.05) +
    geom_point(size = 3) +
    coord_flip() +
    scale_x_discrete(limits = up_order_term) + 
    scale_colour_gradient(low = "red", high = "blue")
}










