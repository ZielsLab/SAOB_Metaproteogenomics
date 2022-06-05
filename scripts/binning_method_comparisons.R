library(tidyverse)

#################################################################
# Comparing sets of bins from NP manual binning and automatic binning across time-series
#################################################################

# checkM stats for both sets of bins
all_TS_bins_info <- read.table("results/re_polished_binning/binning_comparisons/automatic_binning_checkM_stats.tsv", sep="\t", col.names=c("bin1", "lineage1", "completeness1", "contamination1", "strainHeterogeneity1", "size1", "contigs1", "GC1"))

np_bins_info <- read.table("results/re_polished_binning/binning_comparisons/np_bins_checkm_stats.tsv", sep="\t", col.names=c("bin2", "lineage2", "completeness2", "contamination2", "strainHeterogeneity2", "size2", "contigs2", "GC2"))

# ANI comparisons
ani_comparisons <- read.table("results/re_polished_binning/binning_comparisons/comps.txt", sep="\t", col.names = c("bin1", "bin2", "ANI", "l1", "l2")) %>% 
  select(bin1, bin2, ANI) %>% 
  mutate(bin1 = gsub(".fa", "", bin1)) %>% 
  mutate(bin2 = gsub (".fa", "", bin2))

all_binning_comparisons <- left_join(ani_comparisons, all_TS_bins_info) %>% 
  left_join(np_bins_info) %>% 
  select(-lineage2) %>% 
  arrange(lineage1) %>% 
  filter(ANI > 95) %>% 
  mutate(lineage = lineage1) %>% 
  select(-lineage1)

write.csv(all_binning_comparisons, "results/re_polished_binning/binning_comparisons/all_final_binning_comparisons.csv", quote=FALSE, row.names = FALSE)
