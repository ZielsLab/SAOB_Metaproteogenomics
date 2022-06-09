library(tidverse)
library(readxl)


# Set Paths

  ## path to MetaProSIP peptide files
  pep.path <- "~/Documents/UBC_Research/PostDoc/Elizabeth McDaniel/SAOB_Enrichments/Metaproteomics/MetaProSIP/polyplish_asm/MSP_out/peptides/"

  ## protein sample metadata file 
  meta <- read_xlsx("~/Documents/UBC_Research/PostDoc/Elizabeth McDaniel/SAOB_Enrichments/Metaproteomics/metaproteomics-sample-map.xlsx") %>%
    rename(prep_id = `Prep ID`, sample = `Original ID`) %>%
    mutate(prep_id = as.character(prep_id)) %>%
    select(prep_id, sample, time_hr, isotope) 
 
# Read in peptide fileS and merge with metadata
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

 