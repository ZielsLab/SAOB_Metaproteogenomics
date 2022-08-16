library(mmgenome2)
library(vegan)
library(tidyverse)
library(shiny)
library(Biostrings)

#read in files
path = "~/Documents/UBC_Research/PostDoc/Elizabeth McDaniel/SAOB_Enrichments/Metagenomics/metaFlye_raconx1_medakax1_raconilmx1/" #file path
contigs = "~/Documents/UBC_Research/PostDoc/Elizabeth McDaniel/SAOB_Enrichments/Metagenomics/metaFlye_raconx1_medakax1_raconilmx1/contigs_medaka.consensus_racon.ilm_x1.renamed.fasta" #contigs fasta file
cov <- read.table(file=paste0(path, "SAMPLES_MERGED-COVs.txt"), header = TRUE, sep = "\t") #coverage table from anvio
gene.calls <- read.table(file=paste0(path, "contigs.gene.calls"), header = TRUE, sep = "\t") %>% #taxa hits of scgs from anvio anvi-export-gene-calls
  select(gene_callers_id, contig)
scg <- read.table(file=paste0(path, "contigs.db.hits"), header = TRUE, sep = "\t") %>% #taxa hits of scgs from anvio anvi-run-scg-taxonomy 
  filter(source != "Protista_83") %>%
  select(contig, geneID = gene) 
tax <- read_csv(file=paste0(path, "contigs_megan_tax.csv"))

# read assembly as DNA String, filter short contigs
assembly <- readDNAStringSet(contigs, format = "fasta")
list <- names(assembly)[which(width(assembly) >= 2000)] #filter out contigs less than 2000 bp
assembly <- assembly[list]

# mm object
mm <- mmload(
  assembly = assembly,
  coverage = cov, 
  essential_genes = scg, 
  taxonomy = tax,
  kmer_BH_tSNE = TRUE, 
  kmer_size = 4, 
  perplexity = 40, theta=0.5, num_threads = 4, max_iter = 2000
  )

mmstats(mm)
#saveRDS(mm, file = paste0(path, "mm.RDS"))

# initial plot
mmplot(mm, 
       x = 'tSNE1',
       y = 'tSNE2',
       color_by = "family", 
       color_vector = scales::viridis_pal(option = "plasma")(3),
       #color_scale_log10 = TRUE,
       factor_shape = 'solid', 
       alpha = 0.05,
       #locator = TRUE   #run this with locator on first, to select the points
       ) 



#Bin group 1

  # Initial group selection
    selection1 <- data.frame(tSNE1 = c(24.793, 21.435, 21.83, 30.125, 34.47, 34.47, 29.138, 28.348),
                             tSNE2 = c(-7, -11.684, -15.899, -16.836, -13.089, -10.279, -7.351, -5.946))#selection of data points
    
    mmplot(mm, 
           x = 'tSNE1',
           y = 'tSNE2',
           color_by = "family", 
           color_vector = scales::viridis_pal(option = "plasma")(3),
           #color_scale_log10 = TRUE,
           factor_shape = 'solid', 
           alpha = 0.05,
           #locator = TRUE   #run this with locator on first, to select the points, 
           selection = selection1
    ) 
    

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
             factor_shape = 'solid'
             #color_scale_log10 = TRUE, 
             #x_limits = c(38, 45),
             #y_limits = c(0, 2000),
             locator = TRUE,
             #fixed_size = 5,
             #selection = selection
      )
      
      #subset group further using locator 
      selection1.1 <- data.frame(cov_R2Sept2020 = c(2882.561, 2572.206, 2647.444, 3841.841, 4387.313, 3954.697),
                                 cov_R2Nov2019 = c(1589.723, 1410.795, 1132.079, 1168.209, 1436.602, 1593.164))
      selection1.2 <- data.frame(cov_R2Sept2020 = c(7443.84, 5355.997, 5421.83, 7415.626, 8863.949, 8940.491),
                                 cov_R2Nov2019 = c(233.997, 130.769, 51.627, 51.627, 68.832, 223.674))
      #visualize selections
      mmplot(mm_subset1, 
             x = 'cov_R2Sept2020',
             y = 'cov_R2Nov2019',
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
      mmplot(mm_subset1, 
             x = 'cov_R2Sept2020',
             y = 'cov_R2Nov2019',
             color_by = "species", 
             color_vector = scales::viridis_pal(option = "plasma")(3),
             factor_shape = 'solid',
             #color_scale_log10 = TRUE, 
             #x_limits = c(38, 45),
             #y_limits = c(0, 2000),
             #locator = TRUE,
             #fixed_size = 5,
             selection = selection1.2
      )


      #extract subsets as bins
      mm_subset1.1 <- mmextract(mm_subset1, selection = selection1.1)
      mm_subset1.2 <- mmextract(mm_subset1, selection = selection1.2)

      mmstats(mm_subset1.1)
      mmstats(mm_subset1.2)


      # write.table(mm_subset1 %>% 
      #               select(scaffold) %>% 
      #               mutate(bin = "selection_1"),
      #             file= paste0(path, "selection1.tsv"),
      #             col.names = FALSE,
      #             row.names = FALSE, quote = FALSE, sep = "\t"
      # )
        
        #export bins
        mmexport(mm_subset1.1, assembly = assembly,
                 file = paste0(path, "mmgenome_bins/bin1.1.fa"))
        mmexport(mm_subset1.2, assembly = assembly,
                 file = paste0(path, "mmgenome_bins/bin1.2.fa"))

  #remove bins from assembly to avoid duplicate contig binning
        mm.new2 <- mm %>%
          filter(!(scaffold %in% c(mm_subset1.1$scaffold, mm_subset1.2$scaffold)))

## Selection 2

  #initial plotting
        selection2 <- data.frame(tSNE1 = c(13.498, 7.094, 5.721, 8.695, 16.015, 19.903, 21.047, 24.021, 21.276),
                                 tSNE2 = c(-2.667, -12.621, -18.827, -23.862, -22.34, -12.855, -8.171, -4.775, -1.731))
    mmplot(mm.new2, 
           x = 'tSNE1',
           y = 'tSNE2',
           color_by = "cov_R2Mar2020", 
           color_vector = scales::viridis_pal(option = "plasma")(3),
           color_scale_log10 = TRUE, 
           #locator = TRUE #run with locator to get above selection first
           selection = selection2 
    )


  #extract group
  mm_subset2 <- mmextract(mm.new2, selection = selection2)

  #inspect sub plots
    mmplot_pairs(mm_subset2,
                 variables = c("tSNE1",
                               "tSNE2",
                               "cov_R2Sept2020",
                               "cov_R2Nov2019",
                               "cov_R2Mar2020",
                               "cov_R1Sept2020",
                               "gc"),
                 color_by = "order",
                 alpha = 0.4,
                 size_scale = 0.7,
                 textsize = 4,
                 color_vector = scales::viridis_pal(option = "plasma")(3))


    # chose this plot to inspect
    mmplot(mm_subset2, 
           x = 'cov_R2Sept2020',
           y = 'cov_R2Mar2020',
           color_by = "order", 
           #color_vector = scales::viridis_pal(option = "plasma")(3),
           factor_shape = "solid", 
           color_scale_log10 = TRUE, 
           alpha = 0.5,
           #x_limits = c(38, 45),
           #y_limits = c(0, 2000),
           locator = TRUE
           #fixed_size = 5,
           #selection = selection
    )

  #selection2.1 <- data.frame(cov_R2Sept2020 = c(270.83, 251.49, 317.796, 353.712, 983.619, 1309.624, 1038.016, 651.094, 558.29),
   #            tSNE1 = c(18.706, 14.605, 11.158, 9.377, 10.03, 12.024, 16.664, 16.844, 18.363))
  #selection2.2 <- data.frame(cov_R2Sept2020 = c(-24.785, -40.213, -11.596, 5.526, 178.165, 86.852, 22.647),
  #          tSNE1 = c(19.327, 13.543, 8.788, 6.828, 12.961, 14.656, 19.891))
  
   selection2.1 <- data.frame(cov_R2Sept2020 = c(269.64, 334.29, 505.3, 680.481, 966.192, 1176.826, 1251.903, 1306.125, 1296.884, 890.583, 598.218, 316.602),
                              cov_R2Mar2020 = c(11.828, 106.687, 214.937, 214.937, 217.169, 203.777, 159.138, 113.383, 60.931, 19.64, 0.668, -7.144)) 
    selection2.2 <- data.frame(cov_R2Sept2020 = c(8.955, -20.242, -18.902, -6.176, 0.188, 106.243, 318.354, 216.541),
                               cov_R2Mar2020 = c(654.635, 490.911, 305.234, 235.039, 194.28, 224.849, 331.274, 532.801))
      #visualize selections
  mmplot(mm_subset2, 
         x = 'cov_R2Sept2020',
         y = 'cov_R2Mar2020',
         color_by = "order", 
         #color_vector = scales::viridis_pal(option = "plasma")(3),
         factor_shape = "solid", 
         color_scale_log10 = TRUE, 
         alpha = 0.5,
         #x_limits = c(38, 45),
         #y_limits = c(0, 2000),
         #locator = TRUE
         #fixed_size = 5,
         selection = selection2.1
  )
  mmplot(mm_subset2, 
         x = 'cov_R2Sept2020',
         y = 'cov_R2Mar2020',
         color_by = "order", 
         #color_vector = scales::viridis_pal(option = "plasma")(3),
         factor_shape = "solid", 
         color_scale_log10 = TRUE, 
         alpha = 0.5,
         #x_limits = c(38, 45),
         #y_limits = c(0, 2000),
         #locator = TRUE
         #fixed_size = 5,
         selection = selection2.2
  )
    
    #extract subsets as bins
        mm_subset2.1 <- mmextract(mm_subset2, selection = selection2.1)
        mm_subset2.2 <- mmextract(mm_subset2, selection = selection2.2)
        
        mmstats(mm_subset2.1)
        mmstats(mm_subset2.2)
      
        
        #inspect sub plots
        mmplot_pairs(mm_subset2.2,
                     variables = c("tSNE1",
                                   "tSNE2",
                                   "cov_R2Sept2020",
                                   "cov_R2Nov2019",
                                   "cov_R2Mar2020",
                                   "cov_R1Mar2020",
                                   "cov_R1Sept2020",
                                   "gc"),
                     color_by = "order",
                     alpha = 0.4,
                     size_scale = 0.7,
                     textsize = 4,
                     color_vector = scales::viridis_pal(option = "plasma")(3))
    
        
        
        #export bins
        mmexport(mm_subset2.1, assembly = assembly,
                 file = paste0(path, "mmgenome_bins/bin2.1.fa"))
        mmexport(mm_subset2.2, assembly = assembly,
                 file = paste0(path, "mmgenome_bins/bin2.2.fa"))
        
    #remove bins from assembly to avoid duplicate contig binning
        mm.new3 <- mm.new2 %>%
          filter(!(scaffold %in% c(mm_subset2.1$scaffold, mm_subset2.2$scaffold)))
        
## Selection 3
  #initial plotting
   selection3 <- data.frame(tSNE1 = c(0.572, -0.549, 3.151, 6.291, 5.618),
                                 tSNE2 = c(-20.426, -22.683, -23.846, -22.205, -19.537))
    mmplot(mm.new3, 
               x = 'tSNE1',
               y = 'tSNE2',
               color_by = "cov_R2Sept2020", 
               color_vector = scales::viridis_pal(option = "plasma")(3),
               color_scale_log10 = TRUE, 
               #locator = TRUE #run with locator to get above selection first
               selection = selection3
    )
        
        #extract group
        mm_subset3 <- mmextract(mm.new3, selection = selection3)
        
        #inspect sub plots
        mmplot_pairs(mm_subset3,
                     variables = c("tSNE1",
                                   "tSNE2",
                                   "cov_R2Sept2020",
                                   "cov_R2Nov2019",
                                   "cov_R2Mar2020",
                                   "cov_R1Sept2020",
                                   "gc"),
                     color_by = "gc",
                     alpha = 0.4,
                     size_scale = 0.7,
                     textsize = 4,
                     color_vector = scales::viridis_pal(option = "plasma")(3))
        
        # chose this plot to inspect
        mmplot(mm_subset3, 
               x = 'cov_R2Sept2020',
               y = 'cov_R1Sept2020',
               color_by = "gc", 
               color_vector = scales::viridis_pal(option = "plasma")(3),
               color_scale_log10 = TRUE, 
               alpha = 0.5,
               #x_limits = c(38, 45),
               #y_limits = c(0, 2000),
               locator = TRUE
               #fixed_size = 5,
               #selection = selection
        )
        
        selection3.1 <- data.frame(cov_R2Sept2020 = c(32.359, 27.731, 30.805, 37.234, 38.849, 36.818, 37.194),
                                   cov_R1Sept2020 = c(30.221, 26.818, 23.44, 23.49, 26.768, 29.793, 31.103))
        
        #extract subsets as bins
        mm_subset3.1 <- mmextract(mm_subset3, selection = selection3.1)
        mmstats(mm_subset3.1)
        
        
        #export bins
        mmexport(mm_subset3.1, assembly = assembly,
                 file = paste0(path, "mmgenome_bins/bin3.fa"))
        
        #remove bins from assembly to avoid duplicate contig binning
        mm.new4 <- mm.new3 %>%
          filter(!(scaffold %in% c(mm_subset3.1$scaffold)))
        
        
## Selection 4
        
  #initial plotting
  selection4 <- data.frame(tSNE1 = c(-5.03, -8.232, -5.03, 0.003, 0.918, -2.056),
                           tSNE2 = c(-12.738, -15.431, -17.656, -17.187, -14.377, -13.089))
  mmplot(mm.new4, 
               x = 'tSNE1',
               y = 'tSNE2',
               color_by = "cov_R2Mar2020", 
               color_vector = scales::viridis_pal(option = "plasma")(3),
               color_scale_log10 = TRUE, 
               #locator = TRUE #run with locator to get above selection first
               selection = selection4
    )
        
  #extract group
  mm_subset4 <- mmextract(mm.new4, selection = selection4)
  
  #inspect sub plots
  mmplot_pairs(mm_subset4,
               variables = c("tSNE1",
                             "tSNE2",
                             "cov_R2Sept2020",
                             "cov_R2Nov2019",
                             "cov_R2Mar2020",
                             "cov_R1Sept2020",
                             "gc"),
               color_by = "gc",
               alpha = 0.4,
               size_scale = 0.7,
               textsize = 4,
               color_vector = scales::viridis_pal(option = "plasma")(3))
  
  # chose this plot to inspect
  mmplot(mm_subset4, 
         x = 'cov_R2Sept2020',
         y = 'cov_R1Sept2020',
         color_by = "gc", 
         color_vector = scales::viridis_pal(option = "plasma")(3),
         color_scale_log10 = TRUE, 
         alpha = 0.5,
         #x_limits = c(38, 45),
         #y_limits = c(0, 2000),
         locator = TRUE
         #fixed_size = 5,
         #selection = selection
  )
  
    selection4.1 <- data.frame(cov_R2Sept2020 = c(66.346, 61.325, 63.993, 69.797, 74.19, 73.286),
                               cov_R1Sept2020 = c(18.34, 14.324, 9.594, 9.772, 13.967, 18.786))
    #extract group
    mm_subset4.1 <- mmextract(mm_subset4, selection = selection4.1)
    mmstats(mm_subset4.1)
    
    #export bins
    mmexport(mm_subset4.1, assembly = assembly,
             file = paste0(path, "mmgenome_bins/bin4.fa"))
    
    #remove bins from assembly to avoid duplicate contig binning
    mm.new5 <- mm.new4 %>%
      filter(!(scaffold %in% c(mm_subset4.1$scaffold)))

    
## Selection 5
    #initial plotting
    selection5 <- data.frame(tSNE1 = c(-16.238, -22.643, -24.015, -21.728, -14.866, -10.291, -2.514, -5.258),
                             tSNE2 = c(-2.902, -2.902, -5.361, -9.108, -10.162, -10.63, -7.351, -3.604))
    mmplot(mm.new5, 
           x = 'tSNE1',
           y = 'tSNE2',
           color_by = "cov_R2Mar2020", 
           color_vector = scales::viridis_pal(option = "plasma")(3),
           color_scale_log10 = TRUE, 
           #locator = TRUE #run with locator to get above selection first
           selection = selection5
    )   
    
    #extract group
    mm_subset5 <- mmextract(mm.new5, selection = selection5)
    
    #inspect sub plots
    mmplot_pairs(mm_subset5,
                 variables = c("tSNE1",
                               "tSNE2",
                               "cov_R2Sept2020",
                               "cov_R2Nov2019",
                               "cov_R2Mar2020",
                               "cov_R1Sept2020",
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
           color_by = "gc", 
           color_vector = scales::viridis_pal(option = "plasma")(3),
           color_scale_log10 = TRUE, 
           alpha = 0.5,
           #x_limits = c(38, 45),
           #y_limits = c(0, 2000),
           locator = TRUE
           #fixed_size = 5,
           #selection = selection
    )
    
    selection5.1 <- data.frame(cov_R2Sept2020 = c(111.658, 104.49, 119.139, 133.165, 135.97),
               cov_R2Mar2020 = c(37.311, 25.889, 15.102, 21.13, 41.118))
    selection5.2 <- data.frame(cov_R2Sept2020 = c(73.322, 67.4, 69.893, 74.568, 97.009, 100.126),
                               cov_R2Mar2020 = c(83.633, 66.183, 43.339, 36.042, 63.327, 87.44))
    selection5.3 <- data.frame(cov_R2Sept2020 = c(-0.547, -0.618, 59.934, 62.133),
                              cov_R2Mar2020 = c(56.665, -3.617, -2.129, 67.106))
    
    #extract group
    mm_subset5.1 <- mmextract(mm_subset5, selection = selection5.1)
    mmstats(mm_subset5.1)
    mm_subset5.2 <- mmextract(mm_subset5, selection = selection5.2)
    mmstats(mm_subset5.2)
    mm_subset5.3 <- mmextract(mm_subset5, selection = selection5.3)
    mmstats(mm_subset5.3)
    
    #inspect sub plots
    mmplot_pairs(mm_subset5.3,
                 variables = c("tSNE1",
                               "tSNE2",
                               "cov_R2Sept2020",
                               "cov_R2Nov2019",
                               "cov_R2Mar2020",
                               "cov_R1Sept2020",
                               "gc"),
                 color_by = "gc",
                 alpha = 0.4,
                 size_scale = 0.7,
                 textsize = 4,
                 color_vector = scales::viridis_pal(option = "plasma")(3))
    
    # chose this plot to inspect
    mmplot(mm_subset5.3, 
           x = 'tSNE1',
           y = 'tSNE2',
           color_by = "cov_R1Sept2020", 
           color_vector = scales::viridis_pal(option = "plasma")(3),
           color_scale_log10 = TRUE, 
           alpha = 0.5,
           #x_limits = c(38, 45),
           #y_limits = c(0, 2000),
           locator = TRUE
           #fixed_size = 5,
           #selection = selection
    )
    
    selection5.3.1 <- data.frame(tSNE1 = c(-13.255, -14.44, -13.617, -12.464, -11.279, -10.489, -10.654, -11.312),
               tSNE2 = c(-4.24, -4.593, -5.398, -5.662, -5.464, -5.012, -4.372, -4.119))
    
      #re-extract bin
      mm_subset5.3 <- mmextract(mm_subset5.3, selection = selection5.3.1)
      mmstats(mm_subset5.3)
    
    
  #export bins
  mmexport(mm_subset5.1, assembly = assembly,
             file = paste0(path, "mmgenome_bins/bin5.1.fa"))
  mmexport(mm_subset5.2, assembly = assembly,
           file = paste0(path, "mmgenome_bins/bin5.2.fa"))
  mmexport(mm_subset5.3, assembly = assembly,
           file = paste0(path, "mmgenome_bins/bin5.3.fa"))
    
  #remove bins from assembly to avoid duplicate contig binning
  mm.new6 <- mm.new5 %>%
    filter(!(scaffold %in% c(mm_subset5.1$scaffold, mm_subset5.2$scaffold, mm_subset5.3$scaffold)))
    
    
## Selection 6
  selection6 <- data.frame(tSNE1 = c(-35.049, -36.408, -28.936, -21.69, -23.048, -29.615),
                           tSNE2 = c(-18.358, -25.033, -26.438, -21.871, -19.646, -17.656))
  mmplot(mm.new6, 
         x = 'tSNE1',
         y = 'tSNE2',
         color_by = "cov_R2Mar2020", 
         color_vector = scales::viridis_pal(option = "plasma")(3),
         color_scale_log10 = TRUE, 
         #locator = TRUE #run with locator to get above selection first
         selection = selection6
  )   

  #extract group
  mm_subset6 <- mmextract(mm.new6, selection = selection6)
  
  #inspect sub plots
  mmplot_pairs(mm_subset6,
               variables = c("tSNE1",
                             "tSNE2",
                             "cov_R2Sept2020",
                             "cov_R2Nov2019",
                             "cov_R2Mar2020",
                             "cov_R1Sept2020",
                             "gc"),
               color_by = "gc",
               alpha = 0.4,
               size_scale = 0.7,
               textsize = 4,
               color_vector = scales::viridis_pal(option = "plasma")(3))  
  
  # chose this plot to inspect
  mmplot(mm_subset6, 
         x = 'cov_R1Sept2020',
         y = 'cov_R2Sept2020',
         color_by = "gc", 
         color_vector = scales::viridis_pal(option = "plasma")(3),
         color_scale_log10 = TRUE, 
         alpha = 0.5,
         #x_limits = c(38, 45),
         #y_limits = c(0, 2000),
         locator = TRUE
         #fixed_size = 5,
         #selection = selection
  )
  
  selection6.1 <-data.frame(cov_R1Sept2020 = c(-3.677, -8.938, -11.654, -0.636, 16.026, 11.861),
             cov_R2Sept2020 = c(37.631, 32.9, 23.839, 22.104, 29.944, 37.013))
  
  #extract group
  mm_subset6.1 <- mmextract(mm.new6, selection = selection6.1)
  mmstats(mm_subset6.1)
  
  #export bins
  mmexport(mm_subset6.1, assembly = assembly,
           file = paste0(path, "mmgenome_bins/bin6.fa"))
  
  #remove bins from assembly to avoid duplicate contig binning
  mm.new7 <- mm.new6 %>%
    filter(!(scaffold %in% c(mm_subset6.1$scaffold)))
  
  
## Selection 7
  selection7 <- data.frame(tSNE1 = c(-8.556, -12.179, -9.235, -1.537, 2.313, 1.407, -3.122),
                           tSNE2 = c(17.005, 14.78, 10.799, 10.213, 11.736, 15.717, 17.59))
    
  mmplot(mm.new7, 
         x = 'tSNE1',
         y = 'tSNE2',
         color_by = "cov_R2Sept2020", 
         color_vector = scales::viridis_pal(option = "plasma")(3),
         color_scale_log10 = TRUE, 
         #locator = TRUE #run with locator to get above selection first
         selection = selection7
  )  
    