# Methanothermobacter Pangenomics with Anvio 

This page documents steps for creating a pangeome of Methanothermobacter references and genomes assembled in this study using Anvi'o. The ultimate goal of this analysis is to identify core and accessory sets of groups of proteins to then look at the protein expression dyanmics of these groups of proteins of the Methanothermobacter in the SIP microcosm experiment. 

## Collecting genome references 

First I searched the GTDB for "Methanothermobacter" and downloaded those genomes with `ncbi-genome-download`. Then I used only genomes with above 90% completeness and less than 10% redundancy, for a total of 24 reference genomes, which are listed in the `methanothermobacter_gtdb_metadata.csv` file. From this study I assembled 4 different Methanothermobacter genomes, two of which `METHANO1` and `METHANO2` are highly enriched in lab-scale bioreactors along with a putative syntrophic acetate oxidizing bacteria `DTU030_1`. I only included three of the four Methanothermobacter genomes in this analysis because the _Methanothermobacter thermoautotrophicus_ genome only has 79% completion and is 229 contigs, and therefore isn't that great of a genome for this analysis. Whereas the other 3 Methanothermobacter genomes were assembled with Nanopore data, have decently high completion, and are on 5 or 6 contigs - so pretty good genomes. 

## Annotation of Methanothermobacter references 

For the 3 Methanothermobacter genomes assembled in this study, I used the MetaPathways2 program for functional annotation and exploring pathways - which does this by calling genes with Prodigal and making functional annotations with Prokka. Since I have to keep these locus tags the same to compare to in the proteomics data, I am keeping these annotations and will have to load these into Anvi'o with `external-gene-calls`. However for the other reference genomes, I can load these into Anvi'o and have the pipeline itself perform annotation. 

## Starting the Anvi'o Pangenomics Pipeline 

### 
### 
### 