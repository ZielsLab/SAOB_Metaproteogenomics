library(mmgenome2)
library(vegan)
library(tidyverse)
library(shiny)
library(Biostrings)
library(Rtsne)
library(viridis)

#################################################################
# Prepare files and mmgenome object
#################################################################


# project path and contigs
repolished_path <- "results/polished_binning/"
contigs <- paste0(repolished_path, "final_np_polished_contigs_modified.fasta")

# SCG gene calls and taxonomy from Anvi'o 
gene_calls <- read.table(file=paste0(repolished_path, "scg-gene-calls.txt"), header=TRUE, sep="\t") %>% 
  select(gene_callers_id, contig)

scg <- read.table(file=paste0(repolished_path, "scg-taxonomy.txt"), header=TRUE, sep="\t") %>% 
  select(gene_callers_id, gene_name)
# combine gene_calls and scg to get the contig that the essential gene is on
scg_table <- left_join(scg, gene_calls) %>% 
  select(contig, gene_name)

# taxonomy from MEGAN 
tax <- read.table(file=paste0(repolished_path, "megan_tax_modf.tsv"), sep="\t", col.names = c("contig", "taxonomy"))

taxonomy_table <- tax %>% 
  separate(taxonomy, into = c("database", "kingdom", "phylum", "class", "order", "family", "genus", "species"), sep=";") %>% 
  filter(phylum != ' ')

# read assembly as DNA String, filter short contigs
assembly <- readDNAStringSet(contigs, format = "fasta")
list <- names(assembly)[which(width(assembly) >= 2000)] #filter out contigs less than 2000 bp
assembly <- assembly[list]

# coverage table from jgi_summarize_contigs of all samples mapped to the nanopore assembly and filter based on the contigs that remain in the assembly
coverage_table <- read.table(file=paste0(repolished_path, "polished_nanopore-depth.txt"), header=TRUE, sep="\t")
# reformat the coverage table to select just the depth columns 
coverage_filtered <- coverage_table %>% 
  select(contigName, (ends_with("sorted.bam"))) %>% 
  filter(contigName %in% list)
colnames(coverage_filtered) <- gsub("trim", "", colnames(coverage_filtered))
colnames(coverage_filtered) <- gsub(".vs.polished_nanopore.sorted.bam", "", colnames(coverage_filtered))
colnames(np_reads_depth_filtered) <- c("contigName", "NPreads")
all_coverage_table <- left_join(coverage_filtered, np_reads_depth_filtered)


# load into mmgenomes
mm <- mmload(
  assembly = assembly,
  coverage = all_coverage_table, 
  essential_genes = scg_table,
  taxonomy = taxonomy_table,
  kmer_BH_tSNE = TRUE, 
  kmer_size = 4, 
  perplexity = 40, theta=0.5, max_iter = 2000
)

mmstats(mm)

#################################################################
# Base plot to make selections from
#################################################################

# initial plot
mmplot(mm, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 

#################################################################
# Selecting Groups 
#################################################################

#################################
# Selection 1 
#################################

selection1 <- data.frame(tSNE1 = c(7.628, 19.427, 18.733, 7.628),
                         tSNE2 = c(-18.121, -18.844, -13.607, -13.246))
# extract selection 1
mm_subset1 <- mmextract(mm, selection = selection1) 
mmstats(mm_subset1)

# plot pairings 
mmplot_pairs(mm_subset1,
             variables = c("tSNE1",
                           "tSNE2",
                           "cov_R2Sept2020",
                           "cov_R2Nov2019",
                           "cov_NPreads",
                           "gc"),
             color_by = "gc",
             alpha = 0.4,
             size_scale = 0.7,
             textsize = 4,
             color_vector = scales::viridis_pal(option = "plasma")(3))

# chose this plot to inspect
mmplot(mm_subset1, 
       x = 'cov_R2Sept2020',
       y = 'cov_NPreads',
       color_by = "species", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       factor_shape = 'solid',
       #color_scale_log10 = TRUE, 
       #x_limits = c(38, 45),
       #y_limits = c(0, 2000),
       locator = TRUE,
       #fixed_size = 5,
       #selection = selection
)
# subset group 1 
selection1.1 <- data.frame(cov_R2Sept2020 = c(26.28, 26.985, 38.143, 36.703),
                           cov_NPreads = c(11.29, 6.674, 6.92, 11.482))

#visualize the subsetted selection
mmplot(mm_subset1, 
       x = 'cov_R2Sept2020',
       y = 'cov_NPreads',
       color_by = "species", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       factor_shape = 'solid',
       #color_scale_log10 = TRUE, 
       #x_limits = c(38, 45),
       #y_limits = c(0, 2000),
       #locator = TRUE,
       #fixed_size = 5,
       selection = selection1.1
)

# extract the subsetted selection as a bin 
mm_subset1.1 <- mmextract(mm_subset1, selection = selection1.1)
mmstats(mm_subset1.1)

#export selection 1.1
mmexport(mm_subset1.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin1.1.fa"))
