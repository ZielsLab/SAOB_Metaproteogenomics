library(tidyverse)
library(readxl)
library(PNWColors)
library(MetBrewer)
library(Peptides)
library(ggpubr)

# Set Paths

  ## path to MetaProSIP peptide files
  pep.path <- "raw_data/metaproteomics_results_v2/MSP_out/peptides/"

  ## protein sample metadata file 
  meta <- read_xlsx("raw_data/metaproteomics_results_v2/metaproteomics-sample-map.xlsx") %>%
    dplyr::rename(prep_id = `Prep ID`, sample = `Original ID`) %>%
    mutate(prep_id = as.character(prep_id)) %>%
    select(prep_id, sample, time_hr, isotope) 
 
  
# Read in peptide files and merge with metadata
 pep <- tibble(
    file = list.files(pep.path)
  ) %>%
   mutate(file = paste0(pep.path, file)) %>%
   mutate(peptides = map(file, ~ read_tsv(.x, col_names = TRUE) ))

 pept <- bind_rows(pep$peptides) %>%
   dplyr::rename(peptide = `Peptide Sequence`,
          prep_id = `Sample Name`,
          proteins = `Protein Accessions`,
          global_LR = `Global Peptide LR`) %>%
   mutate(prep_id = sapply(strsplit(prep_id, "CP_P"), `[`, 2)) %>%
   mutate(prep_id = sapply(strsplit(prep_id, "_10Jan22"), `[`, 1)) %>%
   mutate(prep_id = as.character(prep_id)) %>%
   left_join(meta, by = "prep_id") 
 
 pept %>% 
   mutate(mw = mw(peptide, monoisotopic = FALSE, avgScale = "expasy", label = "none", aaShift = NULL)) %>% 
   select(peptide, mw)
  
 # Filter for unique peptides
 pept <- pept %>%
   filter(Unique == 1)

 
# plot total peptide labeling distribution
 pal <- pnw_palette("Sunset",1000)
 
 ggplot(pept, aes(x = prep_id, y = global_LR)) + 
   geom_bin_2d() + 
   scale_fill_gradientn(colours=rev(pal)) + 
   facet_grid(cols = vars(time_hr), scales = "free_x") + 
   theme_bw() + 
   ylab("Peptide Label Ratio (%)") + 
   xlab("Sample") + 
   theme(axis.text.x = element_blank())
 
 
# summarize results at protein level
 prot <- pept %>% 
   group_by(proteins, prep_id, time_hr) %>% 
   summarise(n_peptides = n(),
             mean_lr = mean(global_LR), 
             median_lr = median(global_LR),
             mean_RIA1 = mean(`RIA 1`),
             mean_RIA2 = mean(`RIA 2`))
 
# summarize results at MAG level
# read in file of annotations with bin and locus tag columns 
bin_table <- read_tsv("results/pathways/metapathways_annoation_table.txt") %>% 
  mutate(protein = ORF_ID) %>% 
  mutate(MAG = bin) %>% 
  select(MAG, protein)

bin_table_upper <- bin_table %>% 
  mutate(proteins = toupper(protein))

mag <- left_join(bin_table_upper, pept) %>% 
  select(MAG, prep_id, time_hr, global_LR, RIA1 = `RIA 1`, RIA2 = `RIA 2`) %>% 
  drop_na() %>% 
  group_by(MAG, prep_id, time_hr) %>% 
  summarise(n_peptides = n(), 
            mean_lr = mean(global_LR), 
            median_lr = median(global_LR),
            mean_RIA1 = mean(RIA1), 
            mean_RIA2 = mean(RIA2))

 #plot contig LR
 ggplot(mag, aes(x = prep_id, y = mean_lr)) + 
   geom_bin_2d() + 
   scale_fill_gradientn(colours=rev(pal)) + 
   facet_grid(cols = vars(time_hr), scales = "free_x") + 
   theme_bw()
 
 #Plot by MAG
  # number of labelled peptides
 ggplot(mag, aes(x = prep_id, y = MAG)) + 
   geom_tile(aes(fill = log10(n_peptides))) + 
   scale_fill_gradientn(colours=rev(pal)) + 
   facet_grid(cols = vars(time_hr), scales = "free_x") + 
   theme_bw()
 
  # Mean LR
 ggplot(mag, aes(x = prep_id, y = MAG)) + 
   geom_tile(aes(fill = mean_lr)) + 
   scale_fill_gradientn(colours=rev(pal)) + 
   facet_grid(cols = vars(time_hr), scales = "free_x") + 
   theme_bw()
 
 
####### Total Protein approach 

#read in sample names
 lfq.names <- read_csv(file =  "raw_data/metaproteomics_results_v2/lfq_sample_names.csv") %>%
   mutate(prep_id = sapply(strsplit(file, "_P"), `[`, 2)) %>%
   mutate(prep_id = sapply(strsplit(prep_id, "_10Jan"), `[`, 1)) %>%
   mutate(sample = paste0("abundance_", sample)) %>%
   select(prep_id, sample)
 
 #read in protein quantification file
lfq <- read_csv(file =  "raw_data/metaproteomics_results_v2/Fido_Protein_Quant.csv", 
                 skip = 3) %>%
  filter(n_proteins == 1) %>% #select uniquely mapped proteins
  select(protein, abundance_1:abundance_9) %>%
  gather(key = "sample", value = "lfq", -protein) %>%
  left_join(lfq.names,by = "sample") %>%
  left_join(meta, by = "prep_id")

 #normalize protein abundance, filter proteins observed less than 3 times 
lfq <- lfq %>% 
   left_join( lfq %>% filter(lfq > 0) %>% group_by(protein) %>% summarize(obs_cnt = n()) ) %>% #determine number of occurrences per protein 
   filter(obs_cnt > 2) %>% #filter proteins observed less than 3 times
   left_join( lfq %>% group_by(prep_id) %>% summarize(total_lfq = sum(lfq)) ) %>% #calculate total protein intensity per sample
   mutate(lfq_norm = lfq/total_lfq)  #calculate relative abundance in (g prot_mag/g prot_total)

  # by MAG 
  lfq.mag <- bin_table %>%
    select(MAG, protein) %>% 
    left_join(lfq) %>% 
    group_by(MAG, prep_id, time_hr, sample.y) %>% 
    summarise(n_proteins = n(),
              lfq_norm_cum = sum(lfq_norm)) %>% 
    drop_na()
  
  
    # Mean MAG LFQ
    ggplot(lfq.mag, aes(x = prep_id, y = MAG)) + 
      geom_tile(aes(fill = log10(lfq_norm_cum))) + 
      scale_fill_gradientn(name = "Relative protein ab. (log)", colours=rev(pal)) + 
      facet_grid(cols = vars(time_hr), scales = "free_x") + 
      theme_bw()
    
    # combine with total protein quantification 
    prot.quant <- read_csv(file =  "raw_data/metaproteomics_results_v2/SAOB_SIP_protein_extract_quant.csv") %>%
      left_join(meta, by = "sample") %>%
      right_join(lfq.mag, by = "prep_id") %>%
      mutate(mag.prot = protein_ug / 10 * lfq_norm_cum) %>%  #divide proteins (ug) by 10ml to get mg prot/L
      group_by(MAG, time_hr = time_hr.y) %>%
      summarise(mean_prot_mag = mean(mag.prot),
                stdev_prot_mag = sd(mag.prot) ) %>% #summarize total protein quantification per mag over time
      mutate(MAG = tolower(MAG))
    
    # Total protein quant
    ggplot(prot.quant, aes(x = as.factor(time_hr), y = MAG)) + 
      geom_tile(aes(fill = mean_prot_mag)) + 
      scale_fill_gradientn(name = "Protein (mg/L)", colours=rev(pal)) + 
      theme_bw()
    
    mag.tp <- mag %>% mutate(MAG = tolower(MAG)) %>% #total quant of labeled protein carbon per mag
      drop_na() %>%
      group_by(MAG, time_hr) %>%
      summarise(mean_ria_mag = mean(mean_RIA2)/100,
                stdev_ria_mag = sd(mean_RIA2)/100, 
                mean_lr_mag = mean(mean_lr) ,
                stdev_lr_mag = sd(mean_lr)) %>% #determine mean label ratio and RIA per mag over time
      left_join(prot.quant, by = c("MAG", "time_hr")) %>%
      mutate(mag_lab_prot = mean_prot_mag  * mean_lr_mag * mean_ria_mag) %>%
      mutate(mag_lab_prot_std = mag_lab_prot * sqrt((stdev_ria_mag/3/mean_ria_mag)^2 + (stdev_lr_mag/3/mean_lr_mag)^2 + (stdev_prot_mag/3/mean_prot_mag)^2 ) )
    
    mag.3.tp <- mag.tp %>% filter(MAG %in% c("bin14_1", "bin4_1", "bin4_2"))
    mag.other.tp <- mag.tp %>% filter(!(MAG %in% c("bin14_1", "bin4_1", "bin4_2"))) %>%
      mutate(MAG = "All_others") %>%
      group_by(time_hr) %>%
      summarize(MAG = MAG, time_hr = time_hr, mag_lab_prot = sum(mag_lab_prot))
    
    mag.tp.summ <- rbind(mag.3.tp, mag.other.tp)
    # Plot absolute labeled protein per MAG
    ggplot(mag.tp.summ, aes(x = as.factor(time_hr), y = MAG)) + 
      geom_tile(aes(fill = mag_lab_prot)) + 
      scale_fill_gradientn(name = "Labelled protein conc. (g-13C_prot/L)", colours=rev(pal)) + 
      theme_bw()
    
    pal <- c(met.brewer("Renoir")[12], met.brewer("Renoir")[3], met.brewer("Renoir")[10], met.brewer("Renoir")[7])
    
    
    mag_incorporation_plot <- ggplot( mag.tp.summ %>% rbind(tibble( MAG = unique(mag.tp.summ$MAG), time_hr = 0, mag_lab_prot = 0)),
                                      aes(x = time_hr, y = mag_lab_prot, group = MAG)) +
      geom_line(aes(color = MAG), size = 1) + 
      geom_point(aes(color = MAG)) + 
      geom_ribbon(aes(x = time_hr, ymin = mag_lab_prot - mag_lab_prot_std,
                      ymax = mag_lab_prot + mag_lab_prot_std, group = MAG, fill = MAG), alpha = 0.25) +
      scale_color_manual(values = pal, labels=c("All other groups", "DTU068", "METHANO1", "METHANO2")) + 
      scale_fill_manual(values = pal) + 
      ylab("Labelled protein conc. (mg 13C-prot/L)") + 
      xlab("Time (hr)") + 
      guides(color=guide_legend("Group"), fill = "none") +
      theme_pubr() +
      theme(legend.position = "none", axis.title.x = element_text(face="bold", size=12), axis.title.y=element_text(face="bold", size=12), legend.title = element_text(face="bold", size=12), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10))
    
    mag_incorporation_plot_labels <- mag_incorporation_plot + 
      annotate(geom="text", x=190, y=6, label="Methanothermobacter_1", fontface="bold.italic", size=5, color=c(met.brewer("Renoir")[10])) +
      annotate(geom='text', x=250, y=3.5, label="DTU068_1", fontface="bold.italic", size=5, color=c(met.brewer("Renoir")[3])) +
      annotate(geom="text", x=350, y=2.5, label="Methanothermobacter_2", fontface="bold.italic", size=5, color=c(met.brewer("Renoir")[7])) +
      annotate(geom="text", x=300, y=.7, label="All other genomes", fontface="bold.italic", size=5, color=c(met.brewer("Renoir")[12]))
    
    mag_incorporation_plot_labels
    
# plots with mag names - top mags for relative activity and labelling ratio, and supplementary plots for all MAGs
    
    # find top active mags 
    # more than 50 labelled proteins in the _18 prep id of last timepoint
    top_mags <- lfq.mag %>% 
      filter(prep_id == '_18') %>% 
      select(MAG, n_proteins) %>% 
      group_by(MAG) %>% 
      arrange(desc(n_proteins)) %>% 
      filter(n_proteins > 50) %>% 
      pull(MAG)
    
    # read in table of MAG names and join with the top 10 names 
    mag_names <- read.csv("results/bin_tables/SAOB_final_bins_info_names_table.csv") %>% 
      select(bin, specific_name) %>% 
      mutate(MAG = gsub("\\.", "_", bin)) %>% 
      mutate(MAG = sub("-", "_", MAG)) %>% 
      select(MAG, specific_name)
    
      time_labels <- c("24 hr", "144 hr", "408 hr")
      names(time_labels) <- c("24", "144", "408")
    
      mag_relative_activity_top <- lfq.mag %>% 
      filter(MAG %in% top_mags) %>% 
      left_join(mag_names) %>% 
      ggplot(aes(x=prep_id, y = fct_rev(specific_name))) + 
      geom_tile(aes(fill = log10(lfq_norm_cum))) + 
      scale_fill_gradientn(name = "Log Relative \nProtein Abundance", colours=rev(pal)) + 
      facet_grid(cols = vars(time_hr), scales="free_x", labeller = labeller(time_hr = time_labels)) + 
      theme_pubr() +
      theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(face = "bold"), strip.text = element_text(size=12, face="bold"), legend.position="bottom")
      
      mag_relative_activity_top
      
      mag_relative_activity_all <- lfq.mag %>% 
        left_join(mag_names) %>% 
        ggplot(aes(x=prep_id, y = fct_rev(specific_name))) + 
        geom_tile(aes(fill = log10(lfq_norm_cum))) + 
        scale_fill_gradientn(name = "Log Relative \nProtein Abundance", colours=rev(pal)) + 
        facet_grid(cols = vars(time_hr), scales="free_x", labeller = labeller(time_hr = time_labels)) + 
        theme_pubr() +
        theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(face = "bold"), strip.text = element_text(size=12, face="bold"), legend.position="bottom")
    

      labelling_ratio_all <- mag %>%
        left_join(mag_names) %>% 
        ggplot(aes(x=prep_id, y= fct_rev(specific_name))) +
        geom_tile(aes(fill=mean_lr)) +
        scale_fill_gradientn(name = "Mean Labelling Ratio", colours=rev(pnw_palette("Starfish", 1000))) + 
        facet_grid(cols = vars(time_hr), scales = "free_x", labeller = labeller(time_hr = time_labels)) + 
        theme_pubr() +
        theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(face = "bold"), strip.text = element_text(size=12, face="bold"), legend.position="bottom")
      
      mean_RIA_all <- mag %>% 
        left_join(mag_names) %>% 
        ggplot(aes(x=prep_id, y=fct_rev(specific_name))) +
        geom_tile(aes(fill=mean_RIA2)) +
        scale_fill_gradientn(name = "Mean RIA", colours=rev(pnw_palette("Shuksan", 1000))) + 
        facet_grid(cols = vars(time_hr), scales = "free_x", labeller = labeller(time_hr = time_labels)) + 
        theme_pubr() +
        theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(face = "bold"), strip.text = element_text(size=12, face="bold"), legend.position="bottom")
      
      mag_lr_top <- mag %>%
        filter(MAG %in% top_mags) %>% 
        left_join(mag_names) %>% 
        ggplot(aes(x= prep_id, y= fct_rev(specific_name))) +
        geom_tile(aes(fill = log10(n_peptides))) + 
        scale_fill_gradientn(name = "Number of Proteins \nLabelled (log10)", colours=rev(pnw_palette("Lake", 1000))) + 
        facet_grid(cols = vars(time_hr), scales = "free_x", labeller = labeller(time_hr = time_labels)) + 
        theme_pubr() +
        theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(face = "bold"), strip.text = element_text(size=12, face="bold"), legend.position="bottom")
    
      
      mag_lr_all <- mag %>% 
        left_join(mag_names) %>% 
        ggplot(aes(x= prep_id, y= fct_rev(specific_name))) +
        geom_tile(aes(fill = log10(n_peptides))) + 
        scale_fill_gradientn(name = "Number of Proteins \nLabelled (log10)", colours=rev(pnw_palette("Lake", 1000))) + 
        facet_grid(cols = vars(time_hr), scales = "free_x", labeller = labeller(time_hr = time_labels)) + 
        theme_pubr() +
        theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(face = "bold"), strip.text = element_text(size=12, face="bold"), legend.position="bottom")
        
      
      mag %>%
        left_join(mag_names) %>% 
        ggplot(aes(x= prep_id, y= fct_rev(specific_name))) +
        geom_tile(aes(fill = log10(n_peptides))) + 
        scale_fill_gradientn(colours=rev(pal)) + 
        facet_grid(cols = vars(time_hr), scales = "free_x") + 
        theme_bw()


      # arrange grids 
      top_mag_activity_lr
      top_mag_activity_lr <- ggarrange(mag_relative_activity_top, mag_lr_top, labels = c("A", "B"))
      all_mag_supplement <- ggarrange(mag_relative_activity_all, mag_lr_all, labelling_ratio_all, mean_RIA_all, labels=c("A", "B", "C", "D"), ncol = 2, nrow=2)
      
      all_mag_supplement
      
      # output figures
      
      # relative activity and lr for top mags 
      mag_relative_activity
      ggsave("figures/MAG_relative_activity.png", mag_relative_activity, width=15, height=20, units=c("cm"))
      
      ggsave("figures/Top_MAGs_activity.png", top_mag_activity_lr, width=35, height=20, units=c("cm"))
      
      # supplement grid for all MAGs - relative activity, # proteins labelled, mean LR, mean RIA 
      ggsave("figures/all_mags_microcosm_dynamics_supplement.png", all_mag_supplement, width=40, height=45, units=c("cm"))
      
      ggsave("figures/Top_MAGs_activity.pdf", top_mag_activity_lr, width=35, height=20, units=c("cm"))
      
      # sip mag incorporation
      ggsave("figures/SIP_MAG_incorporation.png", mag_incorporation_plot, width=20, height=15, units=c("cm"))
 
      # sip mag incorporation figure with labels
      ggsave("figures/SIP_MAG_incorporation_labels.png", mag_incorporation_plot_labels, width=20, height=15, units=c("cm"))
      
      ggsave("figures/SIP_MAG_incorporation_labels.pdf", mag_incorporation_plot_labels, width=20, height=15, units=c("cm"))
      