library(readr)
library(ggplot2)

dat <- read_csv("outputs/chlamy_cd_results.csv")

top10Genes <- dat %>%
  arrange(desc(log2FoldChange)) %>% 
  head(n = 10)

write_csv(top10Genes, "outputs/chlamy_cd_top.csv")

bot10Genes <- dat %>%
  arrange(log2FoldChange) %>% 
  head(n = 10)

write_csv(top10Genes, "outputs/chlamy_cd_bot.csv")

joined_top10 <- bind_rows(top10Genes, bot10Genes) %>%
  mutate(up_down = if_else(log2FoldChange > 0, "UP", "DOWN")) %>%
  arrange(log2FoldChange)

# To organize our bar plot, we can also extract the order of our organized genes
order_joined <- joined_top10 %>% 
  select(gene_id) %>% 
  pull()

joined_top10 %>% 
  ggplot(aes(x = gene_id, y = log2FoldChange, fill = up_down)) +
  geom_errorbar(aes(ymin= log2FoldChange - lfcSE, ymax =log2FoldChange + lfcSE), width = 0.4) +
  geom_col(color = "black", width = 0.8) +
  scale_x_discrete(limits = order_joined) +
  coord_flip()+ 
  labs(y = "log2(FC)", x= "Gene ID", fill = "Regulation") +
  geom_hline(yintercept = 0) +
  theme_light() +
  scale_y_continuous(limits = c(-17,17), breaks = seq(-15, 15, 5), expand = c(0,0)) +
  scale_fill_manual(values = c("red", "forestgreen")) +
  theme(
    panel.border = element_rect(color = "black"),
    panel.grid = element_line(colour = "white"),
    axis.text = element_text(colour = "black"))
