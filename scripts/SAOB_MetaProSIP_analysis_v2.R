library(tidyverse)
library(readxl)
library(PNWColors)


# Set Paths

  ## path to MetaProSIP peptide files
  pep.path <- "results/metaproteomics_results/MSP_out/peptides/"

  ## protein sample metadata file 
  meta <- read_xlsx("metadata/metaproteomics-sample-map_modified.xlsx") %>%
    rename(prep_id = `Prep ID`, sample = `Original ID`) %>%
    mutate(prep_id = as.character(prep_id)) %>%
    select(prep_id, sample, time_hr, isotope) 
 
  
# Read in peptide files and merge with metadata
 pep <- tibble(
    file = list.files(pep.path)
  ) %>%
   mutate(file = paste0(pep.path, file)) %>%
   mutate(peptides = map(file, ~ read_tsv(.x) ))

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

 
# plot total peptide labeling distribution
 pal <- pnw_palette("Sunset",1000)
 
 ggplot(pept, aes(x = prep_id, y = global_LR)) + 
   geom_bin_2d() + 
   scale_fill_gradientn(colours=rev(pal)) + 
   facet_grid(cols = vars(time_hr), scales = "free_x") + 
   theme_bw()
 
 
# summarize results at protein level
 prot <- pept %>% 
   group_by(proteins, prep_id, time_hr) %>% 
   summarise(n_peptides = n(),
             mean_lr = mean(global_LR), 
             median_lr = median(global_LR),
             mean_RIA1 = mean(`RIA 1`),
             mean_RIA2 = mean(`RIA 2`))
 
# summarize results at MAG level
 mag <- pept %>% 
   mutate(MAG = sapply(strsplit(proteins, "~"), `[`, 1)) %>%
   group_by(MAG, prep_id, time_hr) %>% 
   summarise(n_peptides = n(),
             mean_lr = mean(global_LR), 
             median_lr = median(global_LR),
             mean_RIA1 = mean(`RIA 1`),
             mean_RIA2 = mean(`RIA 2`))

 #plot contig LR
 ggplot(mag, aes(x = prep_id, y = mean_lr)) + 
   geom_bin_2d() + 
   scale_fill_gradientn(colours=rev(pal)) + 
   facet_grid(cols = vars(time_hr), scales = "free_x") + 
   theme_bw()
 
 #Plot by MAG
  # number of labelled peptides
 ggplot(mag, aes(x = prep_id, y = MAG)) + 
   geom_tile(aes(fill = log(n_peptides))) + 
   scale_fill_gradientn(colours=rev(pal)) + 
   facet_grid(cols = vars(time_hr), scales = "free_x") + 
   theme_bw()
 
  # Mean LR
 ggplot(mag, aes(x = prep_id, y = MAG)) + 
   geom_tile(aes(fill = mean_lr)) + 
   scale_fill_gradientn(colours=rev(pal)) + 
   facet_grid(cols = vars(time_hr), scales = "free_x") + 
   theme_bw()
 

 