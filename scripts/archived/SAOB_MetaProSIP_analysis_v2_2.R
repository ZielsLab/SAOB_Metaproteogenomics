library(tidyverse)
library(readxl)
library(PNWColors)
library(sjmisc)


# Set Paths

  ## path to MetaProSIP peptide files
  pep.path <- "results/archived/metaproteomics_results_v1/MSP_out/peptides/"

  ## protein sample metadata file 
  meta <- read_xlsx("metadata/metaproteomics-sample-map_modified.xlsx") %>%
    rename(prep_id = `Prep ID`, sample = `Original ID`) %>%
    mutate(prep_id = as.character(prep_id)) %>%
    select(prep_id, sample, time_hr, isotope) 
 
  
# Read in peptide files and merge with metadata
 pep <- tibble(
    file = list.files(pep.path)
  ) %>%
   mutate(file = paste0(pep.path, file)) %>%
   mutate(peptides = map(file, ~ read_tsv(.x) ))

 pept <- bind_rows(pep$peptides) %>%
   rename(peptide = `Peptide Sequence`,
          prep_id = `Sample Name`,
          proteins = `Protein Accessions`,
          global_LR = `Global Peptide LR`) %>%
   mutate(prep_id = sapply(strsplit(prep_id, "CP_P"), `[`, 2)) %>%
   mutate(prep_id = sapply(strsplit(prep_id, "_10Jan22"), `[`, 1)) %>%
   mutate(prep_id = as.character(prep_id)) %>%
   left_join(meta, by = "prep_id") 
  
 # Filter for unique peptides
 pept <- pept %>%
   filter(Unique == 1)

 
# plot total peptide labeling distribution
 pal <- pnw_palette("Sunset",1000)
 
 ggplot(pept, aes(x = prep_id, y = global_LR)) + 
   geom_bin_2d() + 
   scale_fill_gradientn(colours=rev(pal)) + 
   facet_grid(cols = vars(time_hr), scales = "free_x") + 
   theme_bw()
 
 
# summarize results at protein level
 prot <- pept %>% 
   group_by(proteins, prep_id, time_hr) %>% 
   summarise(n_peptides = n(),
             mean_lr = mean(global_LR), 
             median_lr = median(global_LR),
             mean_RIA1 = mean(`RIA 1`),
             mean_RIA2 = mean(`RIA 2`))
 
# summarize results at MAG level
 mag <- pept %>%
   mutate(MAG = sapply(strsplit(proteins, "~"), `[`, 1)) %>%
   select(MAG, prep_id, time_hr, global_LR, RIA1 =`RIA 1`, RIA2= `RIA 2`) %>%
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
   geom_tile(aes(fill = log(n_peptides))) + 
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
 lfq.names <- read_csv(file =  "results/metaproteomics_results/lfq_sample_names.csv") %>%
   mutate(prep_id = sapply(strsplit(file, "_P"), `[`, 2)) %>%
   mutate(prep_id = sapply(strsplit(prep_id, "_10Jan"), `[`, 1)) %>%
   mutate(sample = paste0("abundance_", sample)) %>%
   select(prep_id, sample)
 
 #read in protein quantification file
lfq <- read_csv(file =  "results/metaproteomics_results/Fido_Protein_Quant.csv", 
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
  lfq.mag <- lfq %>% #summarize at mag level
    mutate(MAG = sapply(strsplit(protein, "~"), `[`, 1)) %>%
    group_by(MAG, prep_id, time_hr, sample.y) %>% 
    summarise(n_proteins = n(),
              lfq_norm_cum = sum(lfq_norm))
  # join with names 
  colnames(group_names)[1] <- c("MAG")
  lfq.mag.names <- left_join(lfq.mag, group_names) %>% 
    drop_na()
  
    # Mean LR
    mag_tot_prot_activity <- ggplot(lfq.mag.names, aes(x = prep_id, y = fct_rev(specific_name))) + 
      geom_tile(aes(fill = log10(lfq_norm_cum))) + 
      scale_fill_gradientn(name = "Relative protein ab. (log)", colours=rev(pal)) + 
      facet_grid(cols = vars(time_hr), scales="free_x") + 
      theme_bw()
    
    # find top 10 active mags
    top_10_mags <- lfq.mag.names %>% 
      filter(prep_id == '_18') %>% 
      select(MAG, n_proteins) %>% 
      group_by(MAG) %>% 
      arrange(desc(n_proteins)) %>% 
      filter(n_proteins > 50) %>% 
      pull(MAG)
    
    prot_activity_top10 <- lfq.mag.names %>% 
      filter(MAG %in% top_10_mags) %>% 
      ggplot(aes(x=prep_id, y = fct_rev(specific_name))) + 
      geom_tile(aes(fill = log10(lfq_norm_cum))) + 
      scale_fill_gradientn(name = "Log Relative Protein Abundance", colours=rev(pal)) + 
      facet_grid(cols = vars(time_hr), scales="free_x") + 
      theme_bw() +
      theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    prot_activity_top10
    
    
    ggsave("figures/mag_tot_prot_activity.png", mag_tot_prot_activity, width=20, height=15, units=c("cm"))
    ggsave("figures/top10_activity.png", prot_activity_top10, width=20, height=10, units=c("cm"))

# combine with total protein quantification 
    prot.quant <- read_csv(file =  "results/metaproteomics_results/SAOB_SIP_protein_extract_quant.csv") %>%
      left_join(meta, by = "sample") %>%
      right_join(lfq.mag, by = "prep_id") %>%
      mutate(mag.prot = protein_ug / 10 * lfq_norm_cum) %>%  #divide by 10ml to get g prot/L
      group_by(MAG, time_hr = time_hr.y) %>%
        summarise(mean_prot_mag = mean(mag.prot),
                  stdev_prot_mag = std(mag.prot) ) %>% #summarize total protein quantification per mag over time
      mutate(MAG = tolower(MAG))
    
        # Total protein quant
        ggplot(prot.quant, aes(x = as.factor(time_hr), y = MAG)) + 
          geom_tile(aes(fill = mean_prot_mag)) + 
          scale_fill_gradientn(name = "Protein (g/L)", colours=rev(pal)) + 
          theme_bw()
    
    mag.tp <- mag %>% mutate(MAG = tolower(MAG)) %>% #total quant of labeled protein carbon per mag
      drop_na() %>%
      group_by(MAG, time_hr) %>%
      summarise(mean_ria_mag = mean(mean_RIA2),
                stdev_ria_mag = std(mean_RIA2), 
                mean_lr_mag = mean(mean_lr) ,
                stdev_lr_mag = std(mean_lr)) %>% #determine mean label ratio and RIA per mag over time
      left_join(prot.quant, by = c("MAG", "time_hr")) %>%
      mutate(mag_lab_prot = mean_prot_mag  * mean_lr_mag * mean_ria_mag)
    

      # Plot absolute labeled protein per MAG
      ggplot(mag.tp, aes(x = as.factor(time_hr), y = MAG)) + 
        geom_tile(aes(fill = mag_lab_prot)) + 
        scale_fill_gradientn(name = "Labelled protein carbon (g-13C_prot/L)", colours=rev(pal)) + 
        theme_bw()
    
      top_mags_incorp <- ggplot( mag.tp %>% rbind(tibble( MAG = unique(mag.tp$MAG), time_hr = 0, mag_lab_prot = 0)) %>%
                filter(MAG %in% c("bin4.1", "bin4.2", "bin14.1")),
              aes(x = time_hr, y = mag_lab_prot, group = MAG)
      ) +  geom_line(aes(color = MAG), size=1.5) + 
        geom_point(aes(color = MAG), size=1.5) +
        theme_bw() +
        theme(legend.position = c("right"))
    
      top_mags_incorp

      ggsave("figures/top3MAGS_isotope_incorporation.png", top_mags_incorp, width=20, height=13, units=c("cm"))
 