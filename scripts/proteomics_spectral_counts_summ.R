library(tidyverse)
library(readxl)

proteomics_sheet <- "raw_data/proteomics/EMSL_51366_Ziels_AJR29053_SC.xlsx"

spectral_counts <- read_xlsx(path = proteomics_sheet, sheet = "sc_protXtab", col_names = TRUE) %>% 
  select(-`Grand Total`) %>% 
  mutate(contigName = str_extract(Protein, "[^_]*_[^_]*")) %>% 
  pivot_longer(cols=starts_with("51366"), names_to="sample", values_to="spectral_count", values_drop_na=TRUE)

sample_map <- read_xlsx(path = "metadata/metaproteomics-sample-map.xlsx", sheet = "Sheet1", col_names = TRUE)

samples_12C <- sample_map %>% 
  filter(isotope=="12C-acetate") %>% 
  pull(run_name)

spectral_counts_12C_samples <- spectral_counts %>%
  filter(sample %in% samples_12C) %>% 
  filter(!grepl('Duplicate proteins', Description)) %>% 
  left_join(anvio_bins_scaffolds) %>% 
  select(bins, contigName, sample, Protein, spectral_count)

total_counts <- spectral_counts_12C_samples %>% 
  mutate_if(is.double, as.integer) %>% 
  group_by(sample) %>% 
  drop_na() %>% 
  summarise(total = sum(spectral_count))

spectral_relative_abundance <- spectral_counts_12C_samples %>% 
  select(bins, Protein, sample, spectral_count) %>% 
  pivot_wider(names_from=sample, values_from=spectral_count) %>% 
  group_by(bins) %>% 
  drop_na() %>% 
  mutate_if(is.double, as.integer) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  pivot_longer(!bins, names_to="sample", values_to="spectral_count") %>% 
  left_join(total_counts) %>% 
  group_by(bins, sample) %>% 
  mutate(rel_abundance = (spectral_count / total) * 100) %>% 
  left_join(saob_bins_groups) %>% 
  select(bins, classification, group, specificGroup, sample, spectral_count, total, rel_abundance, R2Sept2020.sorted.bam) %>% 
  drop_na() %>%
  mutate(phylum = gsub("d__Bacteria;", "", classification)) %>% 
  mutate(phylum = gsub("d__Archaea;", "", phylum)) %>% 
  mutate(specificName = paste0(specificGroup, "_", bins))

spectral_counts_ab_plot <- spectral_relative_abundance %>% 
  ggplot(aes(x=rel_abundance, y=specificName, fill=group)) + 
  geom_boxplot() +
  theme_pubr() +
  scale_fill_brewer(palette="Set2") + 
  xlab("Relative Abudance of Protein Spectral Counts Across Time-Series") + ylab("Genome Name")


ggsave("figures/spectral_counts_rel_abund_boxplot.png", spectral_counts_ab_plot, width=25, height=20, units=c("cm"))

spectral_relative_abundance %>% 
  group_by(sample) %>% 
  summarise(sum = sum(rel_abundance))

ad_saob_anvio_bins_grid <- ggarrange(r2_relative_abundance_plot, spectral_counts_ab_plot, ncol=2, nrow=1, common.legend = TRUE, legend = "bottom", labels=c("A", "B"), widths=c(2.5,2.3))
ad_saob_anvio_bins_grid

ggsave("figures/IWA_AD_abstract/ad_saob_rel_ab_grid.png", ad_saob_anvio_bins_grid, width=40, height=15, units=c("cm"))
