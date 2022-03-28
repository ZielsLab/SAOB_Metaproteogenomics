library(tidyverse)
library(ggpubr)
library(stringr)
library(viridis)
library(readxl)
library(lubridate)

#######################################
# Sampling dates for DNA  
#######################################

sampling_dates <- read_xlsx(path="metadata/EMSL-LS_DNA_Samples_2021_MWM.xlsx", sheet="Sheet1", col_names = TRUE)

sampling_info <- sampling_dates %>%
  mutate(Date = ymd(Date)) %>% 
  mutate(operation_day = Date - first(Date)) %>% 
  mutate(operation_day = as.factor(operation_day))


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
anvio_bins_scaffolds <- read_tsv("results/anvio_binning/all-scaffolds.tsv", col_names = FALSE)
colnames(anvio_bins_scaffolds) <- c("bins", "contigName")
anvio_bins_scaffolds$bins <- gsub("-contigs", "", anvio_bins_scaffolds$bins)

anvio_bins_depth <- left_join(anvio_bins_scaffolds, depth_file)

anvio_sample_depth <- anvio_bins_depth %>% 
  select(bins, contigName, contigLen, R2Sept2020.sorted.bam) %>% 
  drop_na() %>% 
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

anvio_all_samples_depths_table <- left_join(anvio_bins_stats, anvio_all_sample_depths)

anvio_all_samples_depths_info <- left_join(anvio_all_samples_depths_table, sampling_info)

anvio_bins_all_samples_depth_plot <- anvio_all_samples_depths_info %>% 
  ggplot(aes(x=as_factor(operation_day), y=covg, fill=t_phylum)) +
  geom_bar(stat="identity", color="black", size=0.3, width=0.8) +
  theme_pubr() +
  scale_y_continuous(expand=c(0,0)) + 
  theme(axis.text.x = element_text(angle = 85, vjust=0.5, hjust=0.5)) +
  xlab("Operation Day") + ylab("Total Coverage")
anvio_bins_all_samples_depth_plot

ggsave("figures/anvio-prelim-bins-all-samples-depth.png", anvio_bins_all_samples_depth_plot, width=20, height=15, units=c("cm"))

# Preliminary anvio bins check w/ checkM stats and GTDB classification
# Anvio bins from unpolished Nanopore assembly 

checkm_stats <- read_tsv("results/anvio_binning/checkm_stats.tsv")
gtdb_stats <- read_tsv("results/anvio_binning/all_classf_results.tsv")
colnames(gtdb_stats)[1] <- c("Bin")

saob_anvio_table <- left_join(gtdb_stats, checkm_stats)

saob_bins_modf <- saob_anvio_table %>% 
  mutate(bins = gsub("-contigs", "", Bin)) %>% 
  select(bins, classification, Completeness, Contamination, Size, Contigs, GC)

saob_bins_covg_table <- left_join(saob_bins_modf, anvio_sample_depth)

saob_bins_covg_table %>% 
  ggplot(aes(x=Completeness, y=Contamination)) + geom_point(aes(size=total_covg, color=classification)) + theme_pubr()


# Relative abundance for all R2 samples 

total_r2_covg <- anvio_bins_depth %>% 
  select(bins, R2Dec2019.sorted.bam, R2Feb2020.sorted.bam, R2Jan2020.sorted.bam, R2July2020.sorted.bam, R2Mar2020.sorted.bam, R2Nov2019.sorted.bam, R2Sept2020.sorted.bam) %>% 
  mutate_if(is.double, as.integer) %>% 
  pivot_longer(!bins, names_to="sample", values_to="covg") %>% 
  group_by(sample) %>% 
  drop_na() %>% 
  summarise(total = sum(covg))

r2_relative_abundance <- anvio_bins_depth %>% 
  select(bins, R2Dec2019.sorted.bam, R2Feb2020.sorted.bam, R2Jan2020.sorted.bam, R2July2020.sorted.bam, R2Mar2020.sorted.bam, R2Nov2019.sorted.bam, R2Sept2020.sorted.bam) %>% 
  drop_na() %>% 
  mutate_if(is.double, as.integer) %>% 
  group_by(bins) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  pivot_longer(!bins, names_to="sample", values_to="covg") %>% 
  left_join(total_r2_covg) %>% 
  group_by(bins, sample) %>% 
  mutate(rel_abundance = (covg / total) * 100)

write.csv(r2_relative_abundance, "results/anvio_binning/r2_relative_abundance_anvio_bins.csv", quote = FALSE, row.names = FALSE)

anvio_prelim_bins_table <- r2_relative_abundance %>%
  select(bins, sample, rel_abundance) %>% 
  pivot_wider(names_from=sample, values_from = rel_abundance) %>% 
  left_join(saob_bins_modf) %>% 
  drop_na()

write.csv(anvio_prelim_bins_table, "results/anvio_binning/anvio_bins_prelim_info.csv", quote = FALSE, row.names = FALSE)

saob_bins_groups <- read.csv("results/anvio_binning/anvio_bins_prelim_info_groups.csv")
r2_abund_table <- left_join(r2_relative_abundance, saob_bins_groups) %>% 
  drop_na() %>% 
  mutate(sample = gsub(".sorted.bam", "", sample)) %>% 
  mutate(phylum = gsub(";c__.*", "", classification)) %>% 
  mutate(phylum = gsub("d__Bacteria;", "", phylum)) %>% 
  mutate(phylum = gsub("d__Archaea;", "", phylum)) %>% 
  mutate(phylum = str_extract(phylum, "[^_]*_[^_]*_[^_]*"))
  
r2_abund_table$sample <- factor(r2_abund_table$sample, levels=c("R2Nov2019", "R2Dec2019", "R2Jan2020", "R2Feb2020", "R2Mar2020", "R2July2020", "R2Sept2020"))

r2_abund_table_info <- left_join(r2_abund_table, sampling_info)

r2_relative_abundance_plot <- r2_abund_table_info %>% 
  ggplot(aes(x=as_factor(operation_day), y=rel_abundance, fill=group)) +
  geom_bar(stat="identity", color="black", size=0.3, width=0.8) +
  theme_pubr() +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_brewer(palette="Set2") +
  xlab("Operation Day") + ylab("Relative Abundance")
r2_relative_abundance_plot

ggsave("figures/r2_relative_abundance_anvio_bins.png", r2_relative_abundance_plot, width=30, height=20, units=c("cm"))

r2_abund_table_info$code <- paste0(r2_abund_table_info$specificGroup, "_", r2_abund_table_info$bins)

r2_abund_table_info %>% 
  select(code, operation_day, rel_abundance)
  

high_bins_list <- r2_abund_table_info %>% 
  select(bins, sample, code, operation_day, rel_abundance) %>% 
  filter(rel_abundance > 5) %>% 
  pull(bins) %>% 
  unique()

top_heatmap <- r2_abund_table_info %>% 
  filter(bins %in% high_bins_list) %>% 
  ggplot(aes(x=as_factor(operation_day), y=code, fill=rel_abundance)) + 
  geom_tile(color="white") + 
  scale_fill_viridis(option="mako", begin=0, end=1, limits=c(0,40)) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ylab("MAG") +
  xlab("Operation Day")

ggsave("figures/saob_top_lineages_heatmap.png", top_heatmap, width=16, height=10, units=c("cm"))
