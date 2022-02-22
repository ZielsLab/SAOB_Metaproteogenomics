library(tidyverse)
library(stringr)

# bacteria results
hmm_output <- read_tsv("~/Desktop/nanopore_assembly_proteins_sccs_output.tsv", col_names = FALSE)
colnames(hmm_output) <- c("protein", "marker", "pfam", "e_value", "tc")

hmm_output %>% 
  select(protein, marker) %>% 
  mutate(contig = str_extract(protein, "[^_]*_[^_]*")) %>% 
  select(contig, marker) %>% 
  group_by(contig, marker) %>% 
  summarize(sum=n()) %>% 
  filter(sum > 1) %>% 
  arrange(desc(sum))

full_table <- hmm_output %>% 
  select(protein, marker) %>% 
  mutate(contig = str_extract(protein, "[^_]*_[^_]*")) %>% 
  select(contig, marker) %>% 
  group_by(contig, marker) %>% 
  summarize(sum=n()) %>% 
  pivot_wider(id_cols="contig", names_from="marker", values_from="sum")

write_tsv(multiple_copies_on_contig, "results/multiple_copies_on_contig.tsv")
write_tsv(full_table, "results/nanopore_assembly_proteins_SCCs_results.tsv")

# archaea results
arch_hmm_output <- read_tsv("results/nanopore_assembly_proteins_sccs_archaea_output.tsv", col_names = FALSE)
colnames(arch_hmm_output) <- c("protein", "marker", "pfam", "e_value", "tc")

arch_hmm_multiple_copies <- arch_hmm_output %>% 
  select(protein, marker) %>% 
  mutate(contig = str_extract(protein, "[^_]*_[^_]*")) %>% 
  select(contig, marker) %>% 
  group_by(contig, marker) %>% 
  summarize(sum=n()) %>% 
  filter(sum > 1) %>% 
  arrange(desc(sum))

full_arch_table <- arch_hmm_output %>% 
  select(protein, marker) %>% 
  mutate(contig = str_extract(protein, "[^_]*_[^_]*")) %>% 
  select(contig, marker) %>% 
  group_by(contig, marker) %>% 
  summarize(sum=n()) %>% 
  pivot_wider(id_cols="contig", names_from="marker", values_from="sum")

write_tsv(arch_hmm_multiple_copies, "results/archaea_hmm_multiple_copies.tsv")
write_tsv(full_arch_table, "results/archaea_SCCs_table.tsv")
