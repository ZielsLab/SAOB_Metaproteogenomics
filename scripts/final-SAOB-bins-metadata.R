library(tidyverse)
library(ggpubr)
library(viridis)
library(RColorBrewer)
library(lubridate)
library(readxl)

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
# relative abundance plots
#######################################
group_names <- read.csv("results/SAOB_final_bins_info_names_table.csv") %>% 
  select(bin, classification, group_name, specific_name)

# sampling dates
sampling_dates <- read_xlsx(path="metadata/EMSL-LS_DNA_Samples_2021_MWM.xlsx", sheet="Sheet1", col_names = TRUE)

sampling_info <- sampling_dates %>%
  mutate(Date = ymd(Date)) %>% 
  mutate(operation_day = Date - first(Date)) %>% 
  mutate(operation_day = as.factor(operation_day))

r2_abund <- relative_abundance_stats %>% 
  select(starts_with("R2"), bin) %>% 
  pivot_longer(!bin, names_to="sample", values_to="relative_abundance") %>% 
  left_join(group_names) %>% 
  left_join(sampling_info) %>% 
  select(bin, sample, relative_abundance, group_name, operation_day)


r2_abundance_succession <- r2_abund %>% 
  ggplot(aes(x=as_factor(operation_day), y=relative_abundance, fill=group_name)) +
  scale_fill_brewer(palette = "Spectral") +
  geom_bar(stat="identity", color="black", size=0.3, width=0.8) +
  theme_pubr() +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Operation Day") + ylab("Relative Abundance of MAG") +
  theme(legend.position="right", legend.text=element_text(size=20), axis.title.x = element_text(size=15), axis.title.y=element_text(size=15), axis.text.x = element_text(size=15), axis.text.y=element_text(size=15))

r2_abundance_succession

ggsave("figures/SAOB_R2_abundance_succession.png", r2_abundance_succession, width=35, height=17, units=c("cm"))

sept_date <- r2_abund %>% 
  filter(operation_day == '283') %>% 
  ggplot(aes(x=as_factor(operation_day), y=relative_abundance, fill=group)) +
  scale_fill_brewer(palette = "Spectral") +
  geom_bar(stat="identity", color="black", size=0.3, width=0.8) +
  theme_pubr() +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Operation Day") + ylab("Relative Abundance of MAG") +
  theme(legend.position="none")

ggsave("figures/SAOB_experiment_abundance.png", sept_date, width=10, height=15, units=c("cm"))

r2_abund %>% 
  filter(group == 'METHANO1' | group == 'METHANO2') %>% 
  ggplot(aes(x=operation_day, y=fct_rev(group), fill=relative_abundance))  + 
  geom_tile(color="white") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ylab("MAG") +
  xlab("Operation Day")

r2_abund %>% 
  filter(group == 'SYNTROPH1' | group == 'SYNTROPH2') %>% 
  ggplot(aes(x=operation_day, y=fct_rev(group), fill=relative_abundance)) + 
  geom_tile(color="white") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ylab("MAG") +
  xlab("Operation Day")
