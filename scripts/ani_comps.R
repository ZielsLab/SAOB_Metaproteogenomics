library(tidyverse)
library(viridis)

# ANI Comparisons between lineages for DTU and Methanothermobacter to references in the GTDB - ordered by GTDB species number

# DTU 
dtu_comps <- read.table("results/syntroph_comparisons/saob_dtu_comps.txt", header=FALSE, sep="\t") %>% 
  select(-V4, -V5)
colnames(dtu_comps) <- c("genome1", "genome2", "ANI")
dtu_comps$genome1 <- gsub("dtu/ref_genomes/", "", dtu_comps$genome1)
dtu_comps$genome1 <- gsub(".fna", "", dtu_comps$genome1)
dtu_comps$genome2 <- gsub("np_ilm_bins/", "", dtu_comps$genome2)
dtu_comps$genome2 <- gsub(".fa", "", dtu_comps$genome2)

dtu_comps_order <- c("GCA_001513545.1", "GCA_003445655.1", "GCA_012521605.1", "GCA_012516535.1", "GCA_012840405.1", "GCA_012842275.1")

dtu_comps$genome1 <- factor(dtu_comps$genome1, levels=c(dtu_comps_order))

dtu_comps %>% 
  ggplot(aes(x=genome1, y=as_factor(genome2), fill=ANI)) +
  geom_raster() +
  scale_fill_viridis(limits=c(90, 100), option = "D") +
  scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) + 
  theme(axis.text.x= element_text(angle=85, hjust=1))

# Methano
methano_comps <- read.table("results/methano_pangenomics/saob_methano_comps.txt", header=FALSE, sep="\t") %>% 
  select(-V4, -V5)
colnames(methano_comps) <- c("genome1", "genome2", "ANI")
methano_comps$genome1 <- gsub("methano/all_methano_refs/", "", methano_comps$genome1)
methano_comps$genome1 <- gsub(".fna", "", methano_comps$genome1)
methano_comps$genome2 <- gsub("np_ilm_bins/", "", methano_comps$genome2)
methano_comps$genome2 <- gsub(".fa", "", methano_comps$genome2)

methanoA_comps <- methano_comps %>% 
  filter(genome2 == 'bin4.2' | genome2 == 'bin4.1')

methanoA_comps_order <- c("GCA_001507955.1", "GCA_013330535.1",  "GCA_014361345.1", "GCA_003584625.1", "GCA_012840205.1", "GCA_011370395.1", "GCA_012840175.1",  "GCA_013178175.1", "GCA_014361315.1", "GCA_014361335.1", "GCF_003264935.1")

methanoA_comps$genome1 <- factor(methanoA_comps$genome1, levels=c(methanoA_comps_order))

methanoA_comps %>% 
  ggplot(aes(x=genome1, y=as_factor(genome2), fill=ANI)) +
  geom_raster() + 
  scale_fill_viridis(limits=c(80, 100), option = "A") + 
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.text.x= element_text(angle=85, hjust=1))

methano_fam_comps <- methano_comps %>% 
  filter(genome2 == 'R1Nov2019-bin.14' | genome2 == 'bin5.1_modf')

methano_fam_comps_order <- c("GCF_000145295.1", "GCF_008033705.1",  "GCF_009917665.1", "GCF_000828575.1", "GCA_002356395.1", "GCA_012719835.1", "GCA_014361435.1", "GCF_000008645.1", "GCF_003385755.1", "GCF_014889545.1", "GCA_012521115.1", "GCA_900095815.1", "GCF_009914355.1")

methano_fam_comps$genome1 <- factor(methano_fam_comps$genome1, levels=c(methano_fam_comps_order))


methano_fam_comps %>% 
  ggplot(aes(x=genome1, y=as_factor(genome2), fill=ANI)) + 
  geom_raster() + 
  scale_fill_viridis(limits=c(75, 100), option = "B") + 
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.text.x= element_text(angle=85, hjust=1))
