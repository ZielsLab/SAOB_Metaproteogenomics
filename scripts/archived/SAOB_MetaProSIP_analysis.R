library(tidyverse)
library(readxl)


# Set Paths

  ## path to MetaProSIP peptide files
  pep.path <- "results/preliminary_proteomics_results/MSP_out/peptides/"

  ## protein sample metadata file 
  meta <- read_xlsx("metadata/metaproteomics-sample-map_modified.xlsx") %>%
    rename(prep_id = `Prep ID`, sample = `Original ID`) %>%
    mutate(prep_id = as.character(prep_id)) %>%
    select(prep_id, sample, time_hr, isotope) 
 
# Read in peptide fileS and merge with metadata
 pep <- tibble(
    file = list.files(pep.path)
  ) %>%
   mutate(file = paste0(pep.path, file)) %>%
   mutate(peptides = map(file, ~ read.csv(.x) ))

 pept <- bind_rows(pep$peptides) %>%
   rename(peptide = `Peptide Sequence`,
          prep_id = `Sample Name`,
          proteins = `Protein Accessions`,
          global_LR = `Global Peptide LR`) %>%
   mutate(prep_id = sapply(strsplit(prep_id, "CP_P"), `[`, 2)) %>%
   mutate(prep_id = sapply(strsplit(prep_id, "_10Jan22"), `[`, 1)) %>%
   mutate(prep_id = as.character(prep_id)) %>%
   left_join(meta, by = "prep_id") 
  
 # Filter for unique peptides
 pept <- pept %>%
   filter(Unique == 1)

prot.MAG.map <- read_tsv(file="results/re_polished_binning/nanopore_genes_contigs_bins_tax_info.tsv") 
colnames(prot.MAG.map) <- c("proteins", "contig", "bin", "classification")
prot.MAG.map <- prot.MAG.map %>% 
  mutate(proteins = as.character(proteins))

## merge with peptide table
mag <- pept %>%
  left_join(prot.MAG.map, by = "proteins")

mag_data <- pept %>%
  left_join(prot.MAG.map, by = "proteins") %>% 
  group_by(bin, prep_id, time_hr) %>% 
  summarise(n_peptides = n(),
            mean_lr = mean(global_LR), 
            median_lr = median(global_LR),
            mean_RIA1 = mean(`RIA 1`),
            mean_RIA2 = mean(`RIA 2`))

mag %>% 
  ggplot(aes(x=prep_id, y=global_LR)) +
  geom_boxplot() +
  facet_wrap(~ bin)

mag$prep_id <- factor(mag$prep_id, levels=c("_15", "_16", "_17", "_03", "_08", "_10", "_04", "_11", "_12"))
mag_data$prep_id <- factor(mag_data$prep_id, levels=c("_15", "_16", "_17", "_03", "_08", "_10", "_04", "_11", "_12"))

mag %>% 
  ggplot(aes(x=prep_id, y=global_LR)) + 
  geom_boxplot() +
  facet_grid(cols=vars(time_hr), scales="free_x")

mag %>% 
  ggplot(aes(x=prep_id, y=global_LR)) +
  geom_tile() +
  facet_grid(cols=vars(time_hr), rows=vars(bin), scales="free_x")

mag_modf <- mag %>% 
  select(bin, prep_id, time_hr, global_LR, proteins, contig)

new_labels_mags <- mag_modf %>% 
  group_by(bin, prep_id) %>% 
  count()

unbinned_contigs_high <- mag_modf %>% 
  filter(is.na(bin)) %>% 
  filter(global_LR > 0.9) %>% 
  group_by(contig) %>% 
  count() %>% 
  filter(n > 3) %>% 
  arrange(desc(n))


write_tsv(unbinned_contigs_high, "results/re_polished_binning/unbinned_contigs.txt", col_names = FALSE)
