library(tidyverse)
library(ggpubr)
library(viridis)
library(RColorBrewer)

#######################################
# Final SAOB Bins metadata organization
#######################################

# checkm stats 
final_checkm_stats <- read_tsv("results/re_polished_binning/final_bins/checkm_stats.tsv", col_names = TRUE)
colnames(final_checkm_stats) <- c("bin", "lineage", "completeness", "contamination", "heterogeneity", "size", "contigs", "gc")

# GTDB taxonomy 
gtdb_bacteria <- read_tsv("results/re_polished_binning/final_bins/gtdbtk.bac120.summary.tsv", col_names = TRUE) %>% 
  select(user_genome, classification) %>% 
  mutate(bin = user_genome,) %>% 
  select(bin, classification)

gtdb_archaea <- read_tsv("results/re_polished_binning/final_bins/gtdbtk.ar122.summary.tsv", col_names = TRUE) %>% 
  select(user_genome, classification) %>% 
  mutate(bin = user_genome) %>% 
  select(bin, classification)

final_gtdb_stats <- rbind(gtdb_bacteria, gtdb_archaea)

# stats table 
final_stats_table <- left_join(final_gtdb_stats, final_checkm_stats) %>% 
  select(-lineage) %>% 
  mutate(size = size / 1000000)

# write final stats table
write.csv(final_stats_table, "results/SAOB_final_bins_stats_table.csv", quote=FALSE, row.names = FALSE)

# coverm relative abundance 
relative_abundance <- read_tsv("results/re_polished_binning/final_bins/saob_relative_abundance_modf.txt", col_names = TRUE) %>% 
  filter(Genome != "unmapped") %>% 
  mutate(bin = Genome) %>% 
  select(-Genome)

colnames(relative_abundance) <- gsub("-vs-final-bins.sorted.*", "", colnames(relative_abundance))

relative_abundance_stats <- left_join(final_stats_table, relative_abundance) %>%
  select(bin, classification, completeness, contamination, heterogeneity, size, contigs, gc, R1Nov2019, R1Dec2019, R1Jan2020, R1Feb2020, R1Mar2020, R1July2020, R1Sept2020, R2Nov2019, R2Dec2019, R2Jan2020, R2Feb2020, R2Mar2020, R2July2020, R2Sept2020)

relative_abundance_table <- relative_abundance_stats %>% 
  select(-completeness, -contamination, -heterogeneity, -size, -contigs, -gc)

# write out relative abundance table 
write.csv(relative_abundance_table, "results/SAOB_relative_abundance_table.csv", quote=FALSE, row.names = FALSE)

write.csv(relative_abundance_stats, "results/SAOB_final_bins_info_table.csv", quote=FALSE, row.names = FALSE)
  # this needs to be manually edited with the group info for specific figures

relative_abundance_stats %>%
  filter(completeness > 90 & contamination < 10) %>% 
  select(classification)

relative_abundance_stats %>% 
  filter(contigs < 25) %>% 
  select(classification)


#######################################
# Final SAOB Bins metadata organization
#######################################
