library(mmgenome2)
library(vegan)
library(tidyverse)
library(shiny)
library(Biostrings)
library(Rtsne)
library(viridis)

# project path and contigs
repolished_path <- "results/polished_binning/"
contigs <- paste0(repolished_path, "final_np_polished_contigs_modified.fasta")

# coverage table from jgi_summarize_contigs 
coverage_table <- read.table(file=paste0(repolished_path, "r2-np-assembly-samples-depth.txt"), header=TRUE, sep="\t")
  # split covg table for np and ilm reads 
np_reads_depth <- coverage_table %>% 
  select(contigName, final_np_assembly_vs_np_mapping.sorted.bam)
colnames(np_reads_depth) <- c("contig", "depth")
ilm_reads_depth <- coverage_table %>% 
  select(contigName, R2Sept2020.sorted.bam)
colnames(ilm_reads_depth) <- c("contig", "depth")

# SCG gene calls and taxonomy from Anvi'o 
gene_calls <- read.table(file=paste0(repolished_path, "scg-gene-calls.txt"), header=TRUE, sep="\t") %>% 
  select(gene_callers_id, contig)

scg <- read.table(file=paste0(repolished_path, "scg-taxonomy.txt"), header=TRUE, sep="\t") %>% 
  select(gene_callers_id, gene_name)

# taxonomy from MEGAN 
tax <- read.table(file=paste0(repolished_path, "megan_tax_modf.tsv"), sep="\t", col.names = c("contig", "taxonomy"))

taxonomy_table <- tax %>% 
  separate(taxonomy, into = c("database", "kingdom", "phylum", "class", "order", "family", "genus", "species"), sep=";") %>% 
  filter(phylum != ' ')

# read assembly as DNA String, filter short contigs
assembly <- readDNAStringSet(contigs, format = "fasta")
list <- names(assembly)[which(width(assembly) >= 2000)] #filter out contigs less than 2000 bp
assembly <- assembly[list]
ilm_reads_depth_filtered <- ilm_reads_depth %>% 
  filter(contig %in% list)
colnames(ilm_reads_depth_filtered) <- c("scaffold", "coverage")
np_reads_depth_filtered <- np_reads_depth %>% 
  filter(contig %in% list)
colnames(np_reads_depth_filtered) <- c("scaffold", "coverage")

# load into mmgenomes
mm <- mmload(
  assembly = assembly,
  coverage = ilm_reads_depth_filtered, 
  essential_genes = scg,
  taxonomy = taxonomy_table,
  kmer_BH_tSNE = TRUE, 
  kmer_size = 4, 
  perplexity = 40, theta=0.5, max_iter = 2000
)

mmstats(mm)

# initial plot
mmplot(mm, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       #locator = TRUE   #run this with locator on first, to select the points
) 
