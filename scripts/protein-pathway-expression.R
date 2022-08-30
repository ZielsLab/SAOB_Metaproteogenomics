library(tidyverse)
library(Peptides)

# tab-delimited file of each protein sequence
# file generated with fasta-to-tsv.py script

protein_table <- read_tsv("results/pathways_annotations/functional_tables/SAOB-proteins-table.tsv")

protein_table_mw <- protein_table %>% 
  mutate(mw = mw(AA_sequence, monoisotopic = FALSE, avgScale = "expasy", label = "none", aaShift = NULL))

# merge with lfq from total 


# merge with protein annotation files 

