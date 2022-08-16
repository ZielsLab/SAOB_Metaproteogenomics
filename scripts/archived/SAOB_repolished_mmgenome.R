library(mmgenome2)
library(vegan)
library(tidyverse)
library(shiny)
library(Biostrings)
library(Rtsne)
library(viridis)
library(ggpubr)

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

# saving the coverage table for other uses 
coverage_samples <- coverage_table %>% 
  select(contigName, (ends_with("sorted.bam")))
colnames(coverage_samples)[1] <- c("scaffold")

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
full_mmplot <- mmplot(mm, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       # locator = TRUE   #run this with locator on first, to select the points
) 

ggsave("figures/SAOB_mmplot_binning.png", full_mmplot, width=30, height=20, units=c("cm"))

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
mm.new3 <- mm.new2 %>%
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


selection3 <- data.frame(tSNE1 = c(-38.622, -38.641, -28.914, -30.4),
                         tSNE2 = c(20.827, 10.009, 10.417, 21.236))

# extract subset
mm_subset3 <- mmextract(mm, selection = selection3) 
mmstats(mm_subset3)

# visualize pairings   
mmplot_pairs(mm_subset3,
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
mmplot(mm_subset3, 
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
selection3.1 <- data.frame(cov_R2Sept2020 = c(2.769, 2.429, 8.798, 8.568),
                           cov_NPreads = c(3.838, -0.156, -0.2, 3.56))

mm_subset3.1 <- mmextract(mm_subset3, selection = selection3.1)
mmstats(mm_subset3.1)

#export bins
mmexport(mm_subset3.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin3.1.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new4 <- mm.new3 %>%
  filter(!(scaffold %in% c(mm_subset3.1$scaffold)))


#################################
# Selection 4
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new4, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection4 <- data.frame(tSNE1 = c(-25.255, -24.777, -14.892, -16.645),
                         tSNE2 = c(33.279, 24.201, 23.991, 33.657))

# extract subset
mm_subset4 <- mmextract(mm.new4, selection = selection4) 
mmstats(mm_subset4)

# visualize pairings   
mmplot_pairs(mm_subset4,
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
mmplot(mm_subset4, 
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
selection4.1 <- data.frame(cov_R2Sept2020 = c(-4.809, -6.939, 37.594, 40.558),
                           cov_NPreads = c(8.594, -1.562, -1.899, 8.623))

mm_subset4.1 <- mmextract(mm_subset4, selection = selection4.1)
mmstats(mm_subset4.1)

#export bins
mmexport(mm_subset4.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin4.1.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new5 <- mm.new4 %>%
  filter(!(scaffold %in% c(mm_subset4.1$scaffold)))

#################################
# Selection 5
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new5, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection5 <- data.frame(tSNE1 = c(35.171, 36.447, 45.694, 45.338),
                         tSNE2 = c(15.525, 4.547, 4.743, 16.506))

# extract subset
mm_subset5 <- mmextract(mm.new5, selection = selection5) 
mmstats(mm_subset5)

# visualize pairings   
mmplot_pairs(mm_subset5,
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
mmplot(mm_subset5, 
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
selection5.1 <- data.frame(cov_R2Sept2020 = c(-13.081, -13.426, 27.15, 24.422),
                           cov_NPreads = c(11.55, -5.531, -7.098, 12.761))

mm_subset5.1 <- mmextract(mm_subset5, selection = selection5.1)
mmstats(mm_subset5.1)

#export bins
mmexport(mm_subset5.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin5.1.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new6 <- mm.new5 %>%
  filter(!(scaffold %in% c(mm_subset5.1$scaffold)))

#################################
# Selection 6
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new6, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection6 <- data.frame(tSNE1 = c(-35.932, -35.018, -23.527, -23.838),
                         tSNE2 = c(0.234, -11.332, -10.94, -0.158))

# extract subset
mm_subset6 <- mmextract(mm.new6, selection = selection6) 
mmstats(mm_subset6)

# visualize pairings   
mmplot_pairs(mm_subset6,
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
mmplot(mm_subset6, 
       x = 'cov_R2Sept2020',
       y = 'cov_R2Nov2019',
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
selection6.1 <- data.frame(cov_R2Sept2020 = c(-1.04, -1.139, 11.849, 11.079),
                           cov_R2Nov2019 = c(2.307, -1.644, -2.383, 2.389))

mm_subset6.1 <- mmextract(mm_subset6, selection = selection6.1)
mmstats(mm_subset6.1)

#export bins
mmexport(mm_subset6.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin6.1.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new7 <- mm.new6 %>%
  filter(!(scaffold %in% c(mm_subset6.1$scaffold)))

#################################
# Selection 7
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new7, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection7 <- data.frame(tSNE1 = c(-28.732, -29.036, -23.104, -23.713),
                         tSNE2 = c(-22.898, -33.485, -31.72, -23.487))

# extract subset
mm_subset7 <- mmextract(mm.new7, selection = selection7) 
mmstats(mm_subset7)

# visualize pairings   
mmplot_pairs(mm_subset7,
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
mmplot(mm_subset7, 
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
selection7.1 <- data.frame(cov_R2Sept2020 = c(14.963, 14.963, 26.793, 25.922),
                           cov_R2Mar2020 = c(4.912, 0.823, 1.13, 4.883))

mm_subset7.1 <- mmextract(mm_subset7, selection = selection7.1)
mmstats(mm_subset7.1)

#export bins
mmexport(mm_subset7.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin7.1.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new8 <- mm.new7 %>%
  filter(!(scaffold %in% c(mm_subset7.1$scaffold)))

#################################
# Selection 8
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new8, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection8 <- data.frame(tSNE1 = c(24.049, 24.81, 32.263, 31.503),
                         tSNE2 = c(19.838, 9.644, 9.84, 19.838))

# extract subset
mm_subset8 <- mmextract(mm.new8, selection = selection8) 
mmstats(mm_subset8)

# visualize pairings   
mmplot_pairs(mm_subset8,
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
mmplot(mm_subset8, 
       x = 'cov_R2Sept2020',
       y = 'cov_R2Nov2019',
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
selection8.1 <- data.frame(cov_R2Sept2020 = c(5823.804, 5863.811, 8650.976, 8519.721),
                           cov_R2Nov2019 = c(350.388, -12.28, -2.639, 367.214))

selection8.2 <- data.frame(cov_R2Sept2020 = c(2503.21, 2543.217, 5423.732, 5023.661),
                           cov_R2Nov2019 = c(1586.953, 1023.289, 1040.767, 1530.15))

mm_subset8.1 <- mmextract(mm_subset8, selection = selection8.1)
mmstats(mm_subset8.1)

mm_subset8.2 <- mmextract(mm_subset8, selection= selection8.2)
mmstats(mm_subset8.2)

#export bins
mmexport(mm_subset8.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin8.1.fa"))

mmexport(mm_subset8.2, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin8.2.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new9 <- mm.new8 %>%
  filter(!(scaffold %in% c(mm_subset8.1$scaffold, mm_subset8.2$scaffold)))

#################################
# Selection 9
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new9, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection9 <- data.frame(tSNE1 = c(10.055, 9.751, 18.117, 16.444),
                         tSNE2 = c(-11.724, -20.938, -20.154, -11.136))

# extract subset
mm_subset9 <- mmextract(mm.new9, selection = selection9) 
mmstats(mm_subset9)

# visualize pairings   
mmplot_pairs(mm_subset9,
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
mmplot(mm_subset9, 
       x = 'cov_R2Sept2020',
       y = 'cov_R2Nov2019',
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
selection9.1 <- data.frame(cov_R2Sept2020 = c(16.534, 12.24, 62.055, 48.313),
                           cov_R2Nov2019 = c(14.049, -7.501, -9.038, 15.313))


mm_subset9.1 <- mmextract(mm_subset9, selection = selection9.1)
mmstats(mm_subset9.1)


#export bins
mmexport(mm_subset9.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin9.1.fa"))


#remove bins from assembly to avoid duplicate contig binning
mm.new10 <- mm.new9 %>%
  filter(!(scaffold %in% c(mm_subset9.1$scaffold)))

#################################
# Selection 10
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new10, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection10 <- data.frame(tSNE1 = c(20.247, 21.92, 33.024, 30.742),
                          tSNE2 = c(8.86, -5.255, -3.294, 10.036))

# extract subset
mm_subset10 <- mmextract(mm.new10, selection = selection10) 
mmstats(mm_subset10)

# visualize pairings   
mmplot_pairs(mm_subset10,
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
mmplot(mm_subset10, 
       x = 'cov_R2Mar2020',
       y = 'cov_R2Nov2019',
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
selection10.1 <- data.frame(cov_R2Mar2020 = c(10.043, 8.578, 46.18, 36.413),
                            cov_R2Nov2019 = c(12.007, -6.982, -10.58, 14.445))


mm_subset10.1 <- mmextract(mm_subset10, selection = selection10.1)
mmstats(mm_subset10.1)


#export bins
mmexport(mm_subset10.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin10.1.fa"))


#remove bins from assembly to avoid duplicate contig binning
mm.new11 <- mm.new10 %>%
  filter(!(scaffold %in% c(mm_subset10.1$scaffold)))

#################################
# Selection 11
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new11, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection11 <- data.frame(tSNE1 = c(-2.265, -1.809, 3.515, 2.602),
                          tSNE2 = c(-7.215, -14.469, -14.077, -6.823))
# extract subset
mm_subset11 <- mmextract(mm.new11, selection = selection11) 
mmstats(mm_subset11)

# visualize pairings   
mmplot_pairs(mm_subset11,
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
mmplot(mm_subset11, 
       x = 'cov_R2Mar2020',
       y = 'cov_R2Nov2019',
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
selection11.1 <- data.frame(cov_R2Mar2020 = c(44.306, 44.106, 55.618, 55.318),
                            cov_R2Nov2019 = c(1470.209, 931.955, 931.955, 1484.757))


mm_subset11.1 <- mmextract(mm_subset11, selection = selection11.1)
mmstats(mm_subset11.1)


#export bins
mmexport(mm_subset11.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin11.1.fa"))


#remove bins from assembly to avoid duplicate contig binning
mm.new12 <- mm.new11 %>%
  filter(!(scaffold %in% c(mm_subset11.1$scaffold)))

#################################
# Selection 12
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new12, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection12 <- data.frame(tSNE1 = c(10.968, 11.424, 16.9, 16.596),
                          tSNE2 = c(30.228, 21.711, 22.317, 30.591))

# extract subset
mm_subset12 <- mmextract(mm.new12, selection = selection12) 
mmstats(mm_subset12)

# visualize pairings   
mmplot_pairs(mm_subset12,
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
mmplot(mm_subset12, 
       x = 'cov_R2Mar2020',
       y = 'cov_R2Nov2019',
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
selection12.1 <- data.frame(cov_R2Mar2020 = c(5.369, 6.771, 26.751, 23.947),
                            cov_R2Nov2019 = c(88.095, 63.279, 63.87, 89.868))


mm_subset12.1 <- mmextract(mm_subset12, selection = selection12.1)
mmstats(mm_subset12.1)


#export bins
mmexport(mm_subset12.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin12.1.fa"))


#remove bins from assembly to avoid duplicate contig binning
mm.new13 <- mm.new12 %>%
  filter(!(scaffold %in% c(mm_subset12.1$scaffold)))

#################################
# Selection 13
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new13, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection13 <- data.frame(tSNE1 = c(10.055, 9.751, 17.205, 15.531),
                          tSNE2 = c(2.587, -14.273, -12.704, 2.783))

# extract subset
mm_subset13 <- mmextract(mm.new13, selection = selection13) 
mmstats(mm_subset13)

# visualize pairings   
mmplot_pairs(mm_subset13,
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
mmplot(mm_subset13, 
       x = 'cov_NPreads',
       y = 'cov_R2Sept2020',
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
selection13.1 <- data.frame(cov_NPreads = c(77.444, 76.236, 137.859, 139.526),
                            cov_R2Sept2020 = c(1193.13, 467.593, 496.717, 1221.604))

selection13.2 <- data.frame(cov_NPreads = c(-0.901, -2.875, 72.782, 71.173),
                            cov_R2Sept2020 = c(577.378, -35.28, -5.773, 575.393))


mm_subset13.1 <- mmextract(mm_subset13, selection = selection13.1)
mmstats(mm_subset13.1)

mm_subset13.2 <- mmextract(mm_subset13, selection = selection13.2)
mmstats(mm_subset13.2)


#export bins
mmexport(mm_subset13.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin13.1.fa"))

mmexport(mm_subset13.2, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin13.2.fa"))


#remove bins from assembly to avoid duplicate contig binning
mm.new14 <- mm.new13 %>%
  filter(!(scaffold %in% c(mm_subset13.1$scaffold, mm_subset13.2$scaffold)))

#################################
# Selection 14
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new14, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection14 <- data.frame(tSNE1 = c(-2.113, -1.505, 5.644, 4.58),
                          tSNE2 = c(14.741, 8.076, 7.88, 14.741))

# extract subset
mm_subset14 <- mmextract(mm.new14, selection = selection14) 
mmstats(mm_subset14)

# visualize pairings   
mmplot_pairs(mm_subset14,
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
mmplot(mm_subset14, 
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
selection14.1 <- data.frame(cov_R2Sept2020 = c(-2.551, -0.218, 42.91, 33.69),
                            cov_R2Mar2020 = c(48.653, -7.99, -11.251, 46.444))


mm_subset14.1 <- mmextract(mm_subset14, selection = selection14.1)
mmstats(mm_subset14.1)

# further inspect this subsection

mmplot(mm_subset14.1, 
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


selection14.1.1 <- data.frame(cov_R2Sept2020 = c(16.237, 16.374, 26.393, 23.933),
                              cov_NPreads = c(5.976, 1.848, 2.114, 7.185))


mm_subset14.1.1 <- mmextract(mm_subset14, selection = selection14.1.1)
mmstats(mm_subset14.1.1)

#export bins
mmexport(mm_subset14.1.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin14.1.1.fa"))


#remove bins from assembly to avoid duplicate contig binning
mm.new15 <- mm.new14 %>%
  filter(!(scaffold %in% c(mm_subset14.1.1$scaffold)))


#################################
# Selection 15
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new15, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection15 <- data.frame(tSNE1 = c(-14.13, -13.826, -6.676, -8.045),
                          tSNE2 = c(16.113, 8.076, 8.468, 16.506))

# extract subset
mm_subset15 <- mmextract(mm.new15, selection = selection15) 
mmstats(mm_subset15)

# visualize pairings   
mmplot_pairs(mm_subset15,
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
mmplot(mm_subset15, 
       x = 'cov_R2Nov2019',
       y = 'cov_R2Sept2020',
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
selection15.1 <- data.frame(cov_R2Nov2019 = c(-1.341, -1.248, 4.665, 2.201),
                            cov_R2Sept2020 = c(46.367, 24.574, 25.841, 46.874))


mm_subset15.1 <- mmextract(mm_subset15, selection = selection15.1)
mmstats(mm_subset15.1)

#export bins
mmexport(mm_subset15.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin15.1.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new16 <- mm.new15 %>%
  filter(!(scaffold %in% c(mm_subset15.1$scaffold)))


#################################
# Selection 16
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new16, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection16 <- data.frame(tSNE1 = c(-20.823, -19.91, -8.654, -9.11),
                          tSNE2 = c(7.488, 0.626, 0.234, 6.115))

# extract subset
mm_subset16 <- mmextract(mm.new16, selection = selection16) 
mmstats(mm_subset16)

# visualize pairings   
mmplot_pairs(mm_subset16,
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
mmplot(mm_subset16, 
       x = 'cov_R2Nov2019',
       y = 'cov_R2Sept2020',
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
selection16.1 <- data.frame(cov_R2Nov2019 = c(-2.077, -2.281, 2.12, 2.222),
                            cov_R2Sept2020 = c(64.894, -32.325, -43.91, 66.35))

selection16.2 <- data.frame(cov_R2Nov2019 = c(61.901, 63.187, 67.538, 67.276),
                            cov_R2Sept2020 = c(81.097, -18.437, -46.432, 95.947))


mm_subset16.1 <- mmextract(mm_subset16, selection = selection16.1)
mmstats(mm_subset16.1)

mm_subset16.2 <- mmextract(mm_subset16, selection = selection16.2)
mmstats(mm_subset16.2)
  # 45kb, no essential genes - plasmid to save for later??

#export bins
mmexport(mm_subset16.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin16.1.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new17 <- mm.new16 %>%
  filter(!(scaffold %in% c(mm_subset16.1$scaffold)))


#################################
# Selection 17
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new17, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection17 <- data.frame(tSNE1 = c(-18.237, -17.476, -11.24, -12),
                          tSNE2 = c(-2.118, -11.136, -10.548, -1.922))

# extract subset
mm_subset17 <- mmextract(mm.new17, selection = selection17) 
mmstats(mm_subset17)

# visualize pairings   
mmplot_pairs(mm_subset17,
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
mmplot(mm_subset17, 
       x = 'cov_R2Mar2020',
       y = 'cov_R2Sept2020',
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
selection17.1 <- data.frame(cov_R2Mar2020 = c(22.468, 22.399, 29.871, 29.117),
                            cov_R2Sept2020 = c(124.998, 94.842, 96.533, 126.958))


mm_subset17.1 <- mmextract(mm_subset17, selection = selection17.1)
mmstats(mm_subset17.1)

# 45kb, no essential genes - plasmid to save for later??

#export bins
mmexport(mm_subset17.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin17.1.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new18 <- mm.new17 %>%
  filter(!(scaffold %in% c(mm_subset17.1$scaffold)))


#################################
# Selection 18
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new18, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection18 <- data.frame(tSNE1 = c(-7.437, -7.437, -4.091, -4.243),
                          tSNE2 = c(-5.255, -11.136, -10.744, -5.451))

# extract subset
mm_subset18 <- mmextract(mm.new18, selection = selection18) 
mmstats(mm_subset18)

# visualize pairings   
mmplot_pairs(mm_subset18,
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
mmplot(mm_subset18, 
       x = 'cov_R2Mar2020',
       y = 'cov_R2Sept2020',
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
selection18.1 <- data.frame(cov_R2Mar2020 = c(45.444, 46.093, 65.48, 66.556),
                            cov_R2Sept2020 = c(91.553, 59.182, 63.359, 96.774))


mm_subset18.1 <- mmextract(mm_subset18, selection = selection18.1)
mmstats(mm_subset18.1)

#export bins
mmexport(mm_subset18.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin18.1.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new19 <- mm.new18 %>%
  filter(!(scaffold %in% c(mm_subset18.1$scaffold)))


#################################
# Selection 19
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new19, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection19 <- data.frame(tSNE1 = c(-13.065, -13.673, -6.22, -6.829),
                          tSNE2 = c(-0.746, -11.528, -12.116, -1.138))

# extract subset
mm_subset19 <- mmextract(mm.new19, selection = selection19) 
mmstats(mm_subset19)

# visualize pairings   
mmplot_pairs(mm_subset19,
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
mmplot(mm_subset19, 
       x = 'cov_R2Mar2020',
       y = 'cov_R2Sept2020',
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
selection19.1 <- data.frame(cov_R2Mar2020 = c(15.727, 16.147, 22.862, 22.526),
                            cov_R2Sept2020 = c(51.769, 31.917, 33.548, 55.304))


mm_subset19.1 <- mmextract(mm_subset19, selection = selection19.1)
mmstats(mm_subset19.1)

#export bins
mmexport(mm_subset19.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin19.1.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new20 <- mm.new19 %>%
  filter(!(scaffold %in% c(mm_subset19.1$scaffold)))


#################################
# Selection 20
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new20, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection20 <- data.frame(tSNE1 = c(-10.936, -10.327, -5.46, -5.916),
                          tSNE2 = c(-22.898, -30.152, -28.976, -23.094))

# extract subset
mm_subset20 <- mmextract(mm.new20, selection = selection20) 
mmstats(mm_subset20)
  #495 kbp with some essential genes, 29 scaffolds - larger plasmid?

# visualize pairings   
mmplot_pairs(mm_subset20,
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
mmplot(mm_subset20, 
       x = 'cov_R2Mar2020',
       y = 'cov_R2Nov2019',
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
selection20.1 <- data.frame(cov_R2Mar2020 = c(0.213, 0.211, 2.345, 2.308),
                            cov_R2Nov2019 = c(1.476, -1.156, -1.98, 1.96))


mm_subset20.1 <- mmextract(mm_subset20, selection = selection20.1)
mmstats(mm_subset20.1)

#export bins
mmexport(mm_subset20.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin20.1.plasmid.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new21 <- mm.new20 %>%
  filter(!(scaffold %in% c(mm_subset20.1$scaffold)))


#################################
# Selection 21
#################################

# full assembly after 1st bin scaffolds removed
mmplot(mm.new21, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection21 <- data.frame(tSNE1 = c(-2.57, -2.265, 4.123, 2.146),
                          tSNE2 = c(-2.902, -9.372, -8.98, -2.902))

# extract subset
mm_subset21 <- mmextract(mm.new21, selection = selection21) 
mmstats(mm_subset21)

# visualize pairings   
mmplot_pairs(mm_subset21,
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
mmplot(mm_subset21, 
       x = 'cov_R2Mar2020',
       y = 'cov_R2Nov2019',
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
selection21.1 <- data.frame(cov_R2Mar2020 = c(11.949, 14.231, 24.308, 22.787),
                            cov_R2Nov2019 = c(292.011, 212.948, 209.76, 296.633))


mm_subset21.1 <- mmextract(mm_subset21, selection = selection21.1)
mmstats(mm_subset21.1)

#export bins
mmexport(mm_subset21.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin21.1.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new22 <- mm.new21 %>%
  filter(!(scaffold %in% c(mm_subset21.1$scaffold)))


#################################
# Selection 22
#################################

# full assembly after scaffolds removed
mmplot(mm.new22, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 


selection22 <- data.frame(tSNE1 = c(-2.874, -2.417, 7.165, 6.557),
                          tSNE2 = c(15.133, 7.488, 6.704, 15.329))

# extract subset
mm_subset22 <- mmextract(mm.new22, selection = selection22) 
mmstats(mm_subset22)
  # under 1 Mb and doesn't have many essential genes

# stop extracting here after 21 selections to assess quality and 23 total "bins" to assess the quality with CheckM
