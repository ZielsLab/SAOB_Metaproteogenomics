library(mmgenome2)
library(vegan)
library(tidyverse)
library(shiny)
library(Biostrings)
library(Rtsne)
library(viridis)
library(ggpubr)

#################################################################
# Final re polished assembly
# Raconx3 + Medakax3 + Racon Illumina Polishing + Polypolish
# Prepare files and mmgenome object
#################################################################

#project path and contigs
repolished_path <- "results/re_polished_binning/"
contigs <- paste0(repolished_path, "raconx3_medakax3_raconilmx1_polypolish_contigs.fasta")

# SCG HMM hits and gene calls from Anvi'o
gene_calls <- read.table(file=paste0(repolished_path, "saob_gene_calls.txt"), header=TRUE, sep="\t") %>% 
  select(gene_callers_id, contig)

hmm_hits_bacteria <- read.table(file=paste0(repolished_path, "hmm-hits-bacteria.txt"), header=TRUE, sep="\t") %>% 
  select(gene_callers_id, gene_name)

hmm_hits_archaea <- read.table(file=paste0(repolished_path, "hmm-hits-archaea.txt"), header=TRUE, sep="\t") %>% 
  select(gene_callers_id, gene_name)

all_hmm_hits <- rbind(hmm_hits_bacteria, hmm_hits_archaea)

scg_table <- left_join(all_hmm_hits, gene_calls) %>% 
  select(contig, gene_name)

# taxonomy from MEGAN 
tax <- read.table(file=paste0(repolished_path, "saob_contigs_megan_tax_modf.tsv"), sep="\t", col.names = c("contig", "taxonomy"))

taxonomy_table <- tax %>% 
  separate(taxonomy, into = c("database", "kingdom", "phylum", "class", "order", "family", "genus", "species"), sep=";") %>% 
  filter(phylum != ' ')

# read assembly as DNA String, filter short contigs
assembly <- readDNAStringSet(contigs, format = "fasta")
list <- names(assembly)[which(width(assembly) >= 2000)] #filter out contigs less than 2000 bp
assembly <- assembly[list]

# coverage table from jgi_summarize_contigs of all samples mapped to the nanopore assembly and filter based on the contigs that remain in the assembly
coverage_table <- read.table(file=paste0(repolished_path, "r2-np-polished-assembly-samples-depth.txt"), header=TRUE, sep="\t")
# reformat the coverage table to select just the depth columns 
coverage_filtered <- coverage_table %>% 
  select(contigName, (ends_with("sorted.bam"))) %>% 
  filter(contigName %in% list)

colnames(coverage_filtered) <- gsub(".sorted.bam", "", colnames(coverage_filtered))

colnames(coverage_filtered)[1] <- c("scaffold")

# load into mmgenomes
mm <- mmload(
  assembly = assembly,
  coverage = coverage_filtered, 
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
                      locator = TRUE   #run this with locator on first, to select the points
) 


ggsave("figures/SAOB_mmplot_repolished_binning.png", full_mmplot, width=30, height=20, units=c("cm"))

#################################################################
# Selecting Groups 
#################################################################

#################################
# Selection 1 
#################################

# Initial group selection
selection1 <- data.frame(tSNE1 = c(-5.639, -5.523, 0.842, -0.315),
                         tSNE2 = c(-22.591, -34.226, -34.226, -22.954))

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
# ochrobactrum
selection1.1 <- data.frame(cov_R2Sept2020 = c(26.787, 26.591, 37.958, 36.194),
                           cov_R2Nov2019 = c(2.825, -2.253, -3.027, 3.288))

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

selection2 <- data.frame(tSNE1 = c(-1.241, -1.936, 4.314, 3.735),
                         tSNE2 = c(-36.772, -48.589, -48.562, -37.167))

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
selection2.1 <- data.frame(cov_R2Sept2020 = c(17.565, 17.565, 31.448, 33.126),
                           cov_R2Mar2020 = c(3.737, 1.111, 1.578, 3.898))

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


selection3 <- data.frame(tSNE1 = c(-34.456, -34.062, -26.299, -26.534),
                         tSNE2 = c(0.543, -11.15, -9.891, 1.622))

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

# UBA3950
#extract subsets as bins
selection3.1 <- data.frame(cov_R2Sept2020 = c(-11.416, -21.537, 31.028, 27.882),
                           cov_R2Mar2020 = c(64.568, -36.234, -49.126, 63.72))

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


selection4 <- data.frame(tSNE1 = c(-24.233, -23.769, -18.781, -19.013),
                         tSNE2 = c(7.738, -1.436, -0.537, 7.738))

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
selection4.1 <- data.frame(cov_R2Sept2020 = c(5682.482, 5758.447, 8835.038, 8662.359),
                           cov_R2Nov2019 = c(354.278, -13.058, 8.41, 343.837))

selection4.2 <- data.frame(cov_R2Sept2020 = c(2479.282, 2580.569, 5277.334, 5074.76),
                           cov_R2Nov2019 = c(1556.783, 1068.157, 1078.554, 1608.765))

mm_subset4.1 <- mmextract(mm_subset4, selection = selection4.1)
mmstats(mm_subset4.1)

mm_subset4.2 <- mmextract(mm_subset4, selection=selection4.2)
mmstats(mm_subset4.2)

#export bins
mmexport(mm_subset4.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin4.1.fa"))

mmexport(mm_subset4.2, assembly=assembly, file=paste0(repolished_path, "mmgenome_bins/bin4.2.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new5 <- mm.new4 %>%
  filter(!(scaffold %in% c(mm_subset4.1$scaffold, mm_subset4.2$scaffold)))


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


selection5 <- data.frame(tSNE1 = c(9.988, 10.22, 13.7, 14.048),
                         tSNE2 = c(11.696, 5.579, 5.939, 12.595))

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
selection5.1 <- data.frame(cov_R2Sept2020 = c(42.403, 42.725, 54.685, 54.216),
                           cov_R2Mar2020 = c(15.009, 12.815, 12.87, 15.497))

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


selection6 <- data.frame(tSNE1 = c(-28.989, -28.293, -20.173, -21.681),
                         tSNE2 = c(40.477, 32.183, 31.632, 42.469))

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
selection6.1 <- data.frame(cov_R2Sept2020 = c(3.762, 4.082, 9.159, 9.017),
                           cov_R2Mar2020 = c(7.42, 2.331, 2.611, 8.493))

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


selection7 <- data.frame(tSNE1 = c(-29.917, -29.337, -24.929, -25.045),
                         tSNE2 = c(20.508, 12.345, 12.877, 21.218))

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
selection7.1 <- data.frame(cov_R2Sept2020 = c(16.353, 16.556, 30.509, 29.498),
                           cov_R2Nov2019 = c(88.261, 68.734, 69.675, 89.673))

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


selection8 <- data.frame(tSNE1 = c(35.624, 35.972, 42.352, 41.913),
                         tSNE2 = c(9.506, 1.52, 1.52, 10.038))

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
selection8.1 <- data.frame(cov_R2Sept2020 = c(3.421, 2.092, 25.253, 21.006),
                           cov_R2Mar2020 = c(5.876, -0.223, -0.125, 5.755))


mm_subset8.1 <- mmextract(mm_subset8, selection = selection8.1)
mmstats(mm_subset8.1)


#export bins
mmexport(mm_subset8.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin8.1.fa"))


#remove bins from assembly to avoid duplicate contig binning
mm.new9 <- mm.new8 %>%
  filter(!(scaffold %in% c(mm_subset8.1$scaffold)))

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


selection9 <- data.frame(tSNE1 = c(16.136, 16.716, 21.008, 20.312),
                         tSNE2 = c(7.731, -0.432, -0.077, 8.618))

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
selection9.1 <- data.frame(cov_R2Mar2020 = c(-9.833, -8.792, 43.338, 21.835),
                           cov_R2Nov2019 = c(5.841, -1.831, -2.12, 6.535))


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


selection10 <- data.frame(tSNE1 = c(0.914, 0.683, 6.098, 5.176),
                          tSNE2 = c(-9.66, -16.758, -16.403, -9.305))

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
selection10.1 <- data.frame(cov_R2Mar2020 = c(22.306, 22.37, 30.253, 28.651),
                            cov_R2Nov2019 = c(33.468, 19.708, 19.791, 34.469))


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


selection11 <- data.frame(tSNE1 = c(-10.723, -11.529, -2.427, -3.003),
                          tSNE2 = c(33.107, 18.556, 19.798, 34.35))
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
       x = 'cov_R2Sept2020',
       y = 'cov_R2Jan2020',
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
selection11.1 <- data.frame(cov_R2Sept2020 = c(2.541, 6.735, 24.348, 20.993),
                            cov_R2Jan2020 = c(7.248, -2.899, -5.286, 7.544))

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


selection12 <- data.frame(tSNE1 = c(-7.497, -6.575, -2.312, -3.925),
                          tSNE2 = c(7.021, -1.674, -1.852, 7.199))

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
selection12.1 <- data.frame(cov_R2Sept2020 = c(11.733, 11.394, 40.214, 33.433),
                            cov_R2Nov2019 = c(22.566, -15.675, -17.28, 28.865))


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


selection13 <- data.frame(tSNE1 = c(-12.681, -12.681, -8.303, -8.879),
                          tSNE2 = c(-6.998, -13.564, -13.386, -6.288))

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
       x = 'cov_R2Nov2019',
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
selection13.1 <- data.frame(cov_R2Nov2019 = c(1021.187, 1028.559, 1465.945, 1409.429),
                            cov_R2Mar2020 = c(56.442, 41.626, 44.365, 57.065))

mm_subset13.1 <- mmextract(mm_subset13, selection = selection13.1)
mmstats(mm_subset13.1)

#export bins
mmexport(mm_subset13.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin13.1.fa"))



#remove bins from assembly to avoid duplicate contig binning
mm.new14 <- mm.new13 %>%
  filter(!(scaffold %in% c(mm_subset13.1$scaffold)))

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


selection14 <- data.frame(tSNE1 = c(-9.455, -9.225, 15.661, 12.55),
                          tSNE2 = c(16.781, 9.151, 10.925, 21.04))

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
       y = 'gc',
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
selection14.1 <- data.frame(cov_R2Sept2020 = c(507.025, 516.388, 959.586, 923.693),
                            gc = c(50.052, 42.395, 42.618, 50.244))

selection14.2 <- data.frame(cov_R2Sept2020 = c(-23.564, -17.477, 172.138, 191.259),
                            gc = c(52.318, 39.364, 39.651, 52.19))



mm_subset14.1 <- mmextract(mm_subset14, selection = selection14.1)
mmstats(mm_subset14.1)

mm_subset14.2 <- mmextract(mm_subset14, selection = selection14.2)
mmstats(mm_subset14.2)

#export bins
mmexport(mm_subset14.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin14.1.fa"))

mmexport(mm_subset14.2, assembly = assembly, 
         file=paste0(repolished_path, "mmgenome_bins/bin14.2.fa"))


#remove bins from assembly to avoid duplicate contig binning
mm.new15 <- mm.new14 %>%
  filter(!(scaffold %in% c(mm_subset14.1$scaffold)))

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


selection15 <- data.frame(tSNE1 = c(6.444, 7.481, 20.384, 17.043),
                          tSNE2 = c(4.714, -6.998, -8.773, 3.472))

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
selection15.1 <- data.frame(cov_R2Nov2019 = c(-4.641, -6.026, 5.946, 3.718),
                            cov_R2Sept2020 = c(72.021, -43.413, -43.021, 82.013))


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


selection16 <- data.frame(tSNE1 = c(16.121, 16.928, 22.688, 20.039),
                          tSNE2 = c(-32.552, -38.408, -38.585, -30.955))

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
selection16.1 <- data.frame(cov_R2Mar2020 = c(2.159, 2.102, 3.328, 3.414),
                            cov_R2Sept2020 = c(24.811, 17.956, 18.173, 25.605))


mm_subset16.1 <- mmextract(mm_subset16, selection = selection16.1)
mmstats(mm_subset16.1)


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


selection17 <- data.frame(tSNE1 = c(15.084, 14.739, 29.14, 29.831),
                          tSNE2 = c(-13.386, -21.549, -20.662, -10.37))

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
selection17.1 <- data.frame(cov_R2Mar2020 = c(-1.598, -1.117, 2.957, 1.534),
                            cov_R2Sept2020 = c(7.946, 0.584, 1.009, 9.504))


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


selection18 <- data.frame(tSNE1 = c(19.923, 20.039, 31.329, 25.684),
                          tSNE2 = c(19.621, 11.103, 10.038, 19.976))

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
selection18.1 <- data.frame(cov_R2Nov2019 = c(-0.003, -0.003, 0.104, 0.094),
                            cov_R2Sept2020 = c(25.074, 12.126, 12.396, 25.299))


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


selection19 <- data.frame(tSNE1 = c(-1.045, -0.584, 2.757, 3.679),
                          tSNE2 = c(-6.288, -11.789, -11.789, -6.111))

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
selection19.1 <- data.frame(cov_R2Mar2020 = c(16.568, 16.725, 22.555, 21.688),
                            cov_R2Sept2020 = c(51.318, 34.742, 34.955, 54.718))


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


selection20 <- data.frame(tSNE1 = c(-13.948, -13.718, -9.455, -9.455),
                          tSNE2 = c(-3.094, -8.418, -8.773, -2.207))

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
selection20.1 <- data.frame(cov_R2Mar2020 = c(14.843, 15.286, 22.025, 20.766),
                            cov_R2Nov2019 = c(295.526, 197.122, 197.744, 298.54))


mm_subset20.1 <- mmextract(mm_subset20, selection = selection20.1)
mmstats(mm_subset20.1)

#export bins
mmexport(mm_subset20.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin20.1.fa"))

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


selection21 <- data.frame(tSNE1 = c(-8.073, -8.188, 13.126, 12.665),
                          tSNE2 = c(17.136, 11.458, 10.925, 19.443))

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



#export bins
mmexport(mm_subset21, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin21.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new22 <- mm.new21 %>%
  filter(!(scaffold %in% c(mm_subset21$scaffold)))

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


selection22 <- data.frame(tSNE1 = c(-6.921, -6.229, -1.275, -1.275),
                          tSNE2 = c(-4.159, -12.677, -12.144, -3.981))

# extract subset
mm_subset22 <- mmextract(mm.new22, selection = selection22) 
mmstats(mm_subset22)

# chose this plot to inspect
mmplot(mm_subset22, 
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
selection22.1 <- data.frame(cov_R2Mar2020 = c(44.698, 45.152, 67.319, 69.405),
                            cov_R2Sept2020 = c(96.398, 57.919, 58.743, 98.047))


mm_subset22.1 <- mmextract(mm_subset22, selection = selection22.1)
mmstats(mm_subset22.1)

#export bins
mmexport(mm_subset20.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin22.1.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new23 <- mm.new22 %>%
  filter(!(scaffold %in% c(mm_subset22.1$scaffold)))


#################################
# Selection 23
#################################

# full assembly after scaffolds removed
mmplot(mm.new23, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 

# everything remaining is <300kbp in size, plasmids?

selection23 <- data.frame(tSNE1 = c(10.015, 10.246, 14.739, 14.624),
                          tSNE2 = c(13.41, 5.069, 5.069, 13.232))

# extract subset
mm_subset23 <- mmextract(mm.new23, selection = selection23) 
mmstats(mm_subset23) # 500 kb with ~20 essential genes on it? 

# chose this plot to inspect
mmplot(mm_subset23, 
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
selection23.1 <- data.frame(cov_R2Mar2020 = c(5.798, 5.917, 7.003, 7.068),
                            cov_R2Sept2020 = c(24.421, 52.001, 52.187, 23.746))

mm_subset23.1 <- mmextract(mm_subset23, selection = selection23.1)
mmstats(mm_subset23.1)

#export bins
mmexport(mm_subset23.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin23.1.element.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new24 <- mm.new23 %>%
  filter(!(scaffold %in% c(mm_subset23.1$scaffold)))


#################################
# Selection 23
#################################

# full assembly after scaffolds removed
mmplot(mm.new24, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       locator = TRUE   #run this with locator on first, to select the points
) 

selection24 <- data.frame(tSNE1 = c(-5.308, -5.077, 1.951, 0.683),
                          tSNE2 = c(-22.082, -33.439, -33.972, -21.904))

# extract subset
mm_subset24 <- mmextract(mm.new24, selection = selection24) 
mmstats(mm_subset24) # 600 kb with a few essential genes

# chose this plot to inspect
mmplot(mm_subset24, 
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
selection24.1 <- data.frame(cov_R2Mar2020 = c(-5.271, -3.485, 14.941, 9.72),
                            cov_R2Sept2020 = c(45.681, 29.654, 31.138, 45.681))

mm_subset24.1 <- mmextract(mm_subset24, selection = selection24.1)
mmstats(mm_subset24.1)

#export bins
mmexport(mm_subset24.1, assembly = assembly,
         file = paste0(repolished_path, "mmgenome_bins/bin24.1.element.fa"))

#remove bins from assembly to avoid duplicate contig binning
mm.new25 <- mm.new24 %>%
  filter(!(scaffold %in% c(mm_subset24.1$scaffold)))



#################################
# final remaining 
#################################

# full assembly after scaffolds removed
mmplot(mm.new25, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "class", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       # locator = TRUE   #run this with locator on first, to select the points
) 
