library(tidyverse)
library(ggpubr)
library(viridis)
library(RColorBrewer)


# read in contigs to MAG mapping
contigs_map = anvio_bins_scaffolds
colnames(contigs_map) <- c("MAG", "contig")

head(contigs_map)

# read in proteomics sample map file 
samples <- read_xlsx(path = "metadata/metaproteomics-sample-map_modified.xlsx",
          range = 'A1:E19',
          col_names = TRUE) %>%
  select(sample = `Original ID`, prep = `Prep ID`, time_hr, isotope) %>%
  mutate(prep = paste0("P", prep))

# read in SIPPER results and join with MAG table
pep <- read_xlsx(path = "raw_data/sipper/51366_Ziels_SIPPER_FilterPassingResults.xlsx",
          sheet = "xtab_13CSamples_Labeled", 
          skip = 1,
          col_names = TRUE) %>%
  mutate(contig1 = sapply(strsplit(Protein, "_"), `[`, 1), 
         contig2 = sapply(strsplit(Protein, "_"), `[`, 2), 
         contig = paste0(contig1,"_", contig2)) %>%
  left_join(contigs_map, by = 'contig') %>%
  select(P_15:P_12, Protein, contig, MAG)
head(pep)

# gather long
pep <- gather(pep, key = "prep", value = "AFE", -contig, -MAG, -Protein)
head(pep)

bin_group_names <- saob_bins_groups %>% 
  mutate(MAG = bins) %>% 
  select(MAG, classification, group, specificGroup)

pep_names <- left_join(pep, bin_group_names) %>% 
  mutate(specificName = paste0(specificGroup, "_", MAG))

#summary table
mag_afe <- pep %>% drop_na() %>% group_by(MAG, prep) %>% 
  summarize(ave_AFE = mean(AFE), 
            stdev_AFE = sd(AFE),
            enriched_peptides = n()) %>%
  print(n=100)


# join with sample data and plot

ggplot(full_join(samples, pep %>% drop_na(), by = "prep") %>% filter(isotope == "13C-acetate"),
       aes(x = AFE, y = MAG, group = MAG)) + 
  geom_boxplot(aes(color = MAG)) + 
  geom_point(aes(color = MAG)) + 
  facet_grid(rows = vars(time_hr)) + 
  xlab("13C content of labeled peptides (%)") + 
  theme_bw()

sipper_timeseries_plot <- ggplot(full_join(samples, pep_names %>% drop_na(), by = "prep") %>% filter(isotope == "13C-acetate"),
       aes(x = AFE, y = specificName, group = MAG)) + 
  geom_boxplot(aes(color = group)) + 
  geom_point(aes(color = group)) + 
  facet_grid(rows = vars(time_hr)) + 
  xlab("13C content of labeled peptides (%)") + 
  scale_color_manual(values=c("#66C2A5", "#FC8D62", "#8DA0CB", "#A6D854", "#FFD92F", "#E5C494")) +
  ylab("") + 
  theme_bw() + theme(legend.position = "bottom")
sipper_timeseries_plot

ggsave("figures/sipper_timeseries_plot.png", sipper_timeseries_plot, width=17, height=20, units=c("cm"))
       