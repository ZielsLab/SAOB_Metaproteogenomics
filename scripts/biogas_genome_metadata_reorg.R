library(tidyverse)
library(readxl)

# organizing the biogas genomes metadata 

# accessions from NCBI
bioproject_accession_file <- read.table("results/ref_tree/biogas_accessions.txt", sep = "\t", col.names = c("Assembly", "Isolate", "Taxonomy")) %>% 
  filter(Assembly != "Assembly")

# XLSX metadata from publication 
genome_taxonomy <- read_xlsx("results/ref_tree/biogas_genomes_metadata.xlsx", sheet = "MAGs taxonomy", col_names = TRUE) %>% 
  select(`MAGs ID`, `GTDB-Tk classification`) %>% 
  mutate(Assembly = `MAGs ID`, GTDB = `GTDB-Tk classification`) %>% 
  select(Assembly, GTDB) %>% 
  mutate(Isolate = gsub("METABAT_", "", Assembly)) %>% 
  select(Isolate, GTDB)

genome_qual <- read_xlsx("results/ref_tree/biogas_genomes_metadata.xlsx", sheet = "checkM all MAGs", col_names = TRUE) %>% 
  select(`Bin ID`, `Completeness`, `Contamination`, `Genome size [bp]`, `Number of scaffolds`, `GC`) %>% 
  mutate(Assembly = `Bin ID`, Scaffolds = `Number of scaffolds`, Size = `Genome size [bp]`) %>% 
  mutate(Isolate = gsub("METABAT_", "", Assembly)) %>% 
  select(Isolate, Completeness, Contamination, Size, Scaffolds, GC)

genome_metadata <- left_join(genome_taxonomy, genome_qual)

# combine genbank accessions with metadata 
accessions_metadata <- left_join(bioproject_accession_file, genome_metadata) %>% 
  mutate(NCBI_Taxonomy = Taxonomy, GTDB_Taxonomy = GTDB) %>% 
  select(Assembly, Isolate, NCBI_Taxonomy, GTDB_Taxonomy, Completeness, Contamination, Size, Scaffolds, GC)

# write out to the metadata folder 
write.csv(accessions_metadata, "results/ref_tree/biogas_accessions_metadata.csv", quote = FALSE, row.names = FALSE)
