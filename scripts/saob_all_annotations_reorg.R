library(tidyverse)

# SAOB annotations organization 
# Annotations for ORFs and pathways from Metapathways with multiple databases
# Annotations for ORFs from KOFAM Koala using KEGG HMMs 

# Read in the metapathways files where additional column is name of file to give genome name
metapathways_files_path <- "results/pathways/functional_tables/"
annotation_files <- dir(metapathways_files_path, pattern="*.functional_and_taxonomic_table.txt")
metapathways_annotation_tables <- data_frame(filename = annotation_files) %>%
  mutate(file_contents = map(filename, ~ read_tsv(file.path(metapathways_files_path, .)))
  ) %>%
  unnest()

metapathways_annotations <- metapathways_annotation_tables %>% 
  mutate(bin = gsub(".functional_and_taxonomic_table.txt", "", filename)) %>% 
  select(bin, ORF_ID, Contig_Name, ORF_length, start, end, strand, target, product, taxonomy)

# read in the significant, deduplicated annotations for kofam 
kofam_annotations <- read.table("results/pathways/all-saob-kofam-annotations-sig-modf.txt", sep="\t", col.names = c("ORF_ID", "KO", "KEGG_annotation"))

# join the metapathways and kofam annotations by locus tag, where bin name is already split out from above
all_annotation_table <- left_join(kofam_annotations, metapathways_annotations) %>% 
  select(bin, ORF_ID, Contig_Name, ORF_length, start, end, strand, target, product, KO, KEGG_annotation)

# write the annotation tables

write.table(metapathways_annotations, "results/pathways/metapathways_annoation_table.txt", sep="\t", row.names = FALSE, quote = FALSE)
