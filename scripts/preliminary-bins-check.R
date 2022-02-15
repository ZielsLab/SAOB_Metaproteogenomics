library(tidyverse)
library(ggpubr)

depth_file <- read_tsv("results/preliminary_binning/r2-nanopore-assembly-v1-depth.txt")

bins_contigs_file <- read_tsv("results/preliminary_binning/scaffolds-to-bins.tsv", col_names = FALSE)
colnames(bins_contigs_file) <- c("genomeName", "contigName")

checkm_results <- read_tsv("results/preliminary_binning/all-bins-checkm-parsed.tsv", col_names=FALSE)
colnames(checkm_results) <- c("genomeName", "lineage", "completeness", "contamination", "heterogeneity")

bins_depth_table <- left_join(bins_contigs_file, depth_file)

avg_depth <- bins_depth_table %>% 
  select(genomeName, contigName, contigLen, R2Sept2020.sorted.bam) %>% 
  group_by(genomeName) %>% 
  summarise(totalLen = sum(contigLen), total_covg = sum(R2Sept2020.sorted.bam), rel_abundance = mean(R2Sept2020.sorted.bam)) %>%
  arrange(desc(totalLen))

avg_depth_info <- left_join(checkm_results, avg_depth)

avg_depth_info %>% ggplot(aes(x=completeness, y=contamination)) + 
  geom_point(aes(size=total_covg, color=lineage)) +
  theme_bw()
