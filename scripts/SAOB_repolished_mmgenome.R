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

# SCG HMM hits and gene calls from Anvi'o
gene_calls <- read.table(file=paste0(repolished_path, "scg-gene-calls.txt"), header=TRUE, sep="\t") %>% 
  select(gene_callers_id, contig)

hmm_hits_bacteria <- read.table(file=paste0(repolished_path, "hmm-hits-bacteria.txt"), header=TRUE, sep="\t") %>% 
  select(gene_callers_id, gene_name)

hmm_hits_archaea <- read.table(file=paste0(repolished_path, "hmm-hits-archaea.txt"), header=TRUE, sep="\t") %>% 
  select(gene_callers_id, gene_name)

all_hmm_hits <- rbind(hmm_hits_bacteria, hmm_hits_archaea)

scg_table <- left_join(all_hmm_hits, gene_calls) %>% 
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

# Initial group selection
selection1 <- data.frame(tSNE1 = c(-16.676, -16.453, -7.77, -10.219),
                         tSNE2 = c(-34.187, -46.636, -48.199, -35.345))

# extract subset
mm_subset1 <- mmextract(mm, selection = selection1) 
mmstats(mm_subset1)

# visualize pairings   
mmplot_pairs(mm_subset1,
             variables = c("tSNE1",
                           "tSNE2",
                           "cov_R2Sept2020",
                           "cov_R2Nov2019",
                           "cov_R2Mar2020",
                           "gc"),
             color_by = "gc",
             alpha = 0.4,
             size_scale = 0.7,
             textsize = 4,
             color_vector = scales::viridis_pal(option = "plasma")(3))


# chose this plot to inspect
mmplot(mm_subset1, 
       x = 'cov_R2Sept2020',
       y = 'cov_R2Mar2020',
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

#extract subsets as bins
selection1.1 <- data.frame(cov_R2Sept2020 = c(14.616, 14.202, 39.448, 37.793),
                           cov_R2Mar2020 = c(4.521, 0.715, 0.949, 4.224))

mm_subset1.1 <- mmextract(mm_subset1, selection = selection1.1)
mmstats(mm_subset1.1)

#export bins
mmexport(mm_subset1.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin1.1.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new2 <- mm %>%
  filter(!(scaffold %in% c(mm_subset1.1$scaffold)))

#################################
# Selection 2
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new2, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 

selection2 <- data.frame(tSNE1 = c(-4.43, -3.762, 5.366, 2.695),
                         tSNE2 = c(-28.768, -38.404, -37.934, -28.298))

# extract subset
mm_subset2 <- mmextract(mm, selection = selection2) 
mmstats(mm_subset2)

# visualize pairings   
mmplot_pairs(mm_subset2,
             variables = c("tSNE1",
                           "tSNE2",
                           "cov_R2Sept2020",
                           "cov_R2Nov2019",
                           "cov_R2Mar2020",
                           "gc"),
             color_by = "gc",
             alpha = 0.4,
             size_scale = 0.7,
             textsize = 4,
             color_vector = scales::viridis_pal(option = "plasma")(3))


# chose this plot to inspect
mmplot(mm_subset2, 
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

#extract subsets as bins
selection2.1 <- data.frame(cov_R2Sept2020 = c(24.814, 24.566, 41.94, 41.692),
                           cov_NPreads = c(14.225, 8.691, 8.931, 14.627))

mm_subset2.1 <- mmextract(mm_subset2, selection = selection2.1)
mmstats(mm_subset2.1)

#export bins
mmexport(mm_subset2.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin2.1.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new3 <- mm %>%
  filter(!(scaffold %in% c(mm_subset2.1$scaffold)))

#################################
# Selection 3
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new3, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 
