library(tidyverse)

# gene IDS to contigs 
saob_gene_calls <- read.table("results/re_polished_binning/saob_gene_contigs.txt", header = TRUE)

# contigs to nanopore bins 
saob_contigs_bins <- read.table("results/re_polished_binning/np_bins_to_contigs.txt", header=TRUE)

# nanopore bins taxonomy
saob_np_tax <- read_tsv("results/re_polished_binning/binning_comparisons/all_gtdbtk_np_classf.tsv", col_names = TRUE)

gene_contigs_bins <- left_join(saob_gene_calls, saob_contigs_bins)
gene_contigs_bins$bin <- gsub(".fa", "", gene_contigs_bins$bin)
gene_contigs_bins_tax <- left_join(gene_contigs_bins, saob_np_tax)

write_tsv(gene_contigs_bins_tax, "results/re_polished_binning/nanopore_genes_contigs_bins_tax_info.tsv")
