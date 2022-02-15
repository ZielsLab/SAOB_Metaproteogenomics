#!/bin/bash
#SBATCH --account emsla51366                  
#SBATCH --time 3:0:0 
#SBATCH --nodes 1                      
#SBATCH --ntasks-per-node 16  
#SBATCH --job-name minimap2-mapping 
#SBATCH --error minimap2-mapping.err    
#SBATCH --output minimap2-mapping.out          
#SBATCH --mail-user=eamcdani@mail.ubc.ca   
#SBATCH --mail-type ALL  


# paths 
project_path="/tahoma/emsla51366"
reads_path="${project_path}/polishing/nanopore_lrs_trimmed"
consensus_path="${project_path}/polishing/racon_x1"

singularity exec ../51366-apps.sif minimap2 -t 16 -x map-ont -a ${consensus_path}/assembly_racon_x1.fasta ${reads_path}/AD_SIP_trimmed_combined.fastq > ${consensus_path}/assembly_racon_x1_mapping.sam 