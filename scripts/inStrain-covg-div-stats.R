library(tidyverse)
library(ggpubr)
library(viridis)
library(RColorBrewer)


#######################################
# Final SAOB Bins coverage and intra-population diversity stats
#######################################

# fix this with new mapping files
# Inoculum mapping file! 

# inStrain quick profiles 
# read in all coverage files into one file
coverage_path <- "results/final_bins/quick_profiles/"
files <- dir(coverage_path, pattern="*-genomeCoverage.csv")
saob_coverage_info <- data_frame(filename = files) %>%
  mutate(file_contents = map(filename, ~ read.csv(file.path(coverage_path, .)))
  ) %>%
  unnest()

saob_coverage_table <- saob_coverage_info %>% 
  mutate(sample = gsub("-vs-bins-genomeCoverage.csv", "", filename)) %>% 
  mutate(bin = genome) %>% 
  select(sample, bin, coverage, breadth)

saob_coverage_filtered <- saob_coverage_info %>% 
  filter(breadth > .90 & coverage > 10) %>% 
  mutate(sample = gsub("-vs-bins-genomeCoverage.csv", "", filename)) %>% 
  mutate(bin = genome) %>% 
  select(sample, bin, coverage, breadth)

# queues for full inStrain profile analysis 
inStrain_queues <- saob_coverage_filtered %>% 
  select(sample, bin) %>% 
  mutate(filename = paste(sample, "-vs-bins.sorted.bam", sep="")) %>% 
  mutate(genome = paste(bin, ".fa", sep="")) %>% 
  select(filename, genome)

write_tsv(inStrain_queues, "results/final_bins/inStrain-profiles-queue.txt", col_names = FALSE)


# combine relative abundance and coverage stats 
# relative abundance long format 

# r1 
relative_abundance_r1 <- relative_abundance %>% 
  select(bin, InoculumNov2019, R1Nov2019, R1Dec2019, R1Jan2020, R1Feb2020, R1Mar2020, R1July2020, R1Sept2020) %>% 
  pivot_longer(!bin, names_to="sample", values_to="relative_abundance")
relative_abundance_r1$sample <- factor(relative_abundance_r1$sample, levels=c("InoculumNov2019", "R1Nov2019", "R1Dec2019", "R1Jan2020", "R1Feb2020", "R1Mar2020", "R1July2020", "R1Sept2020"))
r1_abundance_coverage <- left_join(relative_abundance_r1, saob_coverage_table)
r1_abundance_coverage_info <- left_join(r1_abundance_coverage, final_stats_table)
r1_abundance_coverage_info$sample <- factor(r1_abundance_coverage_info$sample, levels=c("InoculumNov2019", "R1Nov2019", "R1Dec2019", "R1Jan2020", "R1Feb2020", "R1Mar2020", "R1July2020", "R1Sept2020"))

# r1 plots
r1_abundance_coverage_info %>% 
  select(bin, classification, sample, relative_abundance) %>% 
  ggplot(aes(x=as_factor(sample), y=relative_abundance, fill=classification)) +
  geom_bar(stat="identity", color="black", size=0.3, width=0.8) +
  theme_pubr() +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Sample") + ylab("Relative Abundance")

r1_abundance_coverage_info %>% 
  select(bin, sample, relative_abundance) %>% 
  ggplot(aes(x=sample, y=bin, fill=relative_abundance))


# r2
relative_abundance_r2 <- relative_abundance %>% 
  select(bin, InoculumNov2019, R2Nov2019, R2Dec2019, R2Jan2020, R2Feb2020, R2Mar2020, R2July2020, R2Sept2020) %>% 
  pivot_longer(!bin, names_to="sample", values_to="relative_abundance")
relative_abundance_r2$sample <- factor(relative_abundance_r2$sample, levels=c("InoculumNov2019", "R2Nov2019", "R2Dec2019", "R2Jan2020", "R2Feb2020", "R2Mar2020", "R2July2020", "R2Sept2020"))
r2_abundance_coverage <- left_join(relative_abundance_r2, saob_coverage_table)
r2_abundance_coverage_info <- left_join(r2_abundance_coverage, final_stats_table)
r2_abundance_coverage_info$sample <- factor(r2_abundance_coverage_info$sample, levels=c("InoculumNov2019", "R2Nov2019", "R2Dec2019", "R2Jan2020", "R2Feb2020", "R2Mar2020", "R2July2020", "R2Sept2020"))

# r2 plots
r2_abundance_coverage_info %>% 
  select(bin, classification, sample, relative_abundance) %>% 
  ggplot(aes(x=as_factor(sample), y=relative_abundance, fill=classification)) +
  geom_bar(stat="identity", color="black", size=0.3, width=0.8) +
  theme_pubr() +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Sample") + ylab("Relative Abundance")

r2_abundance_coverage %>% 
  select(bin, sample, relative_abundance) %>% 
  ggplot(aes(x=as_factor(sample), y=bin, fill=relative_abundance)) +
  geom_tile(color="white") +
  scale_fill_viridis(option="plasma", trans="sqrt") +
  theme_pubr() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  xlab("Sample") + ylab("Genome")
