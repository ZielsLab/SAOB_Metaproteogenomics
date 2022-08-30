library(tidyverse)
library(Peptides)
library(ggpubr)

# tab-delimited file of each protein sequence
# file generated with fasta-to-tsv.py script

protein_table <- read_tsv("results/pathways_annotations/functional_tables/SAOB-proteins-table.tsv")

protein_table_mw <- protein_table %>% 
  mutate(mw = mw(AA_sequence, monoisotopic = FALSE, avgScale = "expasy", label = "none", aaShift = NULL)) %>% 
  mutate(protein = locus_tag) %>% 
  select(protein, AA_sequence, mw)

# merge with lfq from total protein approach in SAOB_MetaProSIP_analysis script
# lfq filtered out proteins observed less than 3 times, lfq_norm (or total protein) is lfq/total_lfq
# protein concentration is lfq_norm/mw for a protein, use this for looking at expression

protein_table_mw_lfq <- left_join(lfq, protein_table_mw) %>% 
  mutate(protein_abundance = lfq_norm / mw)

protein_table_mw_lfq %>%
  arrange(desc(protein_abundance)) %>% 
  head(100) %>% 
  ggplot(aes(x=prep_id, y=protein)) +
  geom_tile(aes(fill=protein_abundance)) + 
  facet_grid(cols = vars(time_hr), scales="free_x")


# merge with protein annotation files 
# metapathways 
metapathways_table <- read_tsv("results/pathways_annotations/metapathways_annoation_table.txt") %>% 
  mutate(protein = ORF_ID) %>% 
  select(bin, protein, ORF_length, product)

protein_table_mw_lfq_metapathways <- left_join(protein_table_mw_lfq, metapathways_table) %>% 
  mutate(protein_tag = paste(product, " (", protein, ")"))

# kegg  
kegg_table <- read_tsv("results/pathways_annotations/all-saob-kofam-annotations-sig-modf.txt", col_names = c("protein", "KO", "annotation"))

protein_table_all_annotations <- left_join(protein_table_mw_lfq_metapathways, kegg_table)

# expression over time-series in top 3 MAGs 
protein_table_mw_lfq_metapathways %>% 
  filter(bin == 'bin14_1' | bin == 'bin4_1' | bin == "bin4_2") %>% 
  group_by(bin) %>% 
  slice_max(order_by = protein_abundance, n=75) %>% 
  ggplot(aes(x=prep_id, y=protein_tag)) +
  geom_tile(aes(fill=protein_abundance)) + 
  facet_grid(cols = vars(time_hr), rows = vars(bin), scales="free") +
  theme_pubr()



# expression of redundant proteins (did not map uniquely)
lfq_redundant <- read_csv(file =  "raw_data/metaproteomics_results_v2/Fido_Protein_Quant.csv", 
                skip = 3) %>%
  filter(n_proteins > 1) %>% #select uniquely mapped proteins
  select(protein, abundance_1:abundance_9) %>%
  gather(key = "sample", value = "lfq", -protein) %>%
  left_join(lfq.names,by = "sample") %>%
  left_join(meta, by = "prep_id")
