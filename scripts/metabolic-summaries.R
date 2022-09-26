library(tidyverse)
library(ggpubr)
library(PNWColors)

### Metabolic summaries of SAOB MAGs 

# Select genomes with relative abundance above 1% in any sample along the time-series 
top_saob_genomes <- read.csv("results/SAOB_final_bins_info_table.csv") %>% 
  select(bin, classification, matches('R2')) %>% 
  filter(if_any(starts_with("R2"), ~. > 1)) %>% 
  pull(bin)

top_saob_genomes_metadata <- read.csv("results/SAOB_final_bins_info_table.csv") %>% 
  select(bin, classification) %>% 
  filter(bin %in% top_saob_genomes)

# all KEGG decoder results for top genomes 
top_kegg_decoder <- read_tsv("results/pathways_annotations/all_kegg_decoder_results_modified.txt") %>% 
  filter(bin %in% top_saob_genomes)

# subset KEGG decoder results for top genomes 
subset_kegg_decoder <- read.csv("results/pathways_annotations/kegg-decoder-subset-pathways.csv") %>% 
  filter(bin %in% top_saob_genomes) %>% 
  left_join(top_saob_genomes_metadata) %>% 
  mutate(bin_name = paste(classification, "(", bin, ")")) %>% 
  select(-bin, -classification) %>% 
  pivot_longer(-bin_name, names_to="pathway", values_to="percent")

subset_kegg_decoder$pathway <- gsub("\\.", " ", subset_kegg_decoder$pathway)

pathway_heatmap <- subset_kegg_decoder %>% 
  ggplot(aes(x=pathway, y=fct_rev(bin_name), fill=percent)) +
  geom_raster() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) + 
  scale_fill_gradientn(colours = pnw_palette("Shuksan2",1000)) +
  theme_pubr() +
  theme(axis.text.y = element_blank(), axis.text.x=element_text(angle = 85, hjust=1), axis.title.y = element_blank(), axis.title.x=element_blank())



ggsave("figures/pathway_heatmap.png", pathway_heatmap, width=15, height=35, units=c("cm"))
