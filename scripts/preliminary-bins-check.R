library(tidyverse)
library(ggpubr)

#######################################
# Preliminary binning stats and coverage 
#######################################

# depth file of contigs
depth_file <- read_tsv("results/preliminary_binning/r2-nanopore-assembly-v1-depth.txt")

# preliminary metabat bins 
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

metabat_bins_plot <- avg_depth_info %>% ggplot(aes(x=completeness, y=contamination)) + 
  geom_point(aes(size=total_covg, color=lineage)) +
  theme_bw() + 
  ggtitle("Preliminary binning with MetaBAT")

#######################################
# Preliminary binning stats and coverage from anvi'o binning
#######################################

anvio_bins_stats <- read_tsv("results/preliminary_binning/anvio_binning/bins_summary.txt")
anvio_bins_scaffolds <- read_tsv("results/preliminary_binning/anvio_binning/scaffolds-to-bins.tsv", col_names = FALSE)
colnames(anvio_bins_scaffolds) <- c("bins", "contigName")
anvio_bins_scaffolds$bins <- gsub("-contigs", "", anvio_bins_scaffolds$bins)

anvio_bins_depth <- left_join(anvio_bins_scaffolds, depth_file)

anvio_sample_depth <- anvio_bins_depth %>% 
  select(bins, contigName, contigLen, R2Sept2020.sorted.bam) %>% 
  group_by(bins) %>% 
  summarise(totalLen = sum(contigLen), total_covg = sum(R2Sept2020.sorted.bam), rel_abundance = mean(R2Sept2020.sorted.bam)) %>% 
  arrange(desc(totalLen))

anvio_depth_info <- left_join(anvio_bins_stats, anvio_sample_depth)

anvio_bins_plot <- anvio_depth_info %>% 
  ggplot(aes(x=percent_completion, y=percent_redundancy)) + 
  geom_point(aes(size=total_covg, color=t_phylum)) +
  theme_bw() + 
  ggtitle("Preliminary binning with Anvi'o")

prelim_bins_qual <- ggarrange(metabat_bins_plot, anvio_bins_plot, ncol=1, labels=c("A", "B"))

ggsave("figures/preliminary-bins-qual-comps.png", prelim_bins_qual, width=20, height=25, units=c("cm"))

anvio_all_sample_depths <- anvio_bins_depth %>% 
  select(bins, ends_with("sorted.bam")) %>% 
  group_by(bins) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  pivot_longer(!bins, names_to="sample", values_to="covg") %>% 
  mutate(sample = gsub(".sorted.bam", "", sample))

anvio_all_samples_depths_info <- left_join(anvio_bins_stats, anvio_all_sample_depths)

anvio_all_samples_depths_info$sample <- factor(anvio_all_samples_depths_info$sample, levels=c("R1Nov2019", "R1Dec2019", "R1Jan2020", "R1Feb2020", "R1Mar2020", "R1July2020", "R1Sept2020", "R2Nov2019", "R2Dec2019", "R2Jan2020", "R2Feb2020", "R2Mar2020", "R2July2020", "R2Sept2020"))

anvio_bins_all_samples_depth_plot <- anvio_all_samples_depths_info %>% 
  ggplot(aes(x=as_factor(sample), y=covg, fill=t_phylum)) +
  geom_bar(stat="identity", color="black", size=0.3, width=0.8) +
  theme_pubr() +
  scale_y_continuous(expand=c(0,0)) + 
  theme(axis.text.x = element_text(angle = 85, vjust=0.5, hjust=0.5)) +
  xlab("Sample") + ylab("Total Coverage")

anvio_bins_all_samples_depth_plot
ggsave("figures/anvio-prelim-bins-all-samples-depth.png", anvio_bins_all_samples_depth_plot, width=20, height=15, units=c("cm"))
