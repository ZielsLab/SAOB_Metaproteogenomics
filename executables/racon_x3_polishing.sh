#!/bin/bash
#SBATCH --account emsla51366                  
#SBATCH --time 5:0:0 
#SBATCH --nodes 2                      
#SBATCH --ntasks-per-node 16  
#SBATCH --job-name racon_x3_polishing
#SBATCH --error racon_x3_polishing.err    
#SBATCH --output racon_x3_polishing.out          
#SBATCH --mail-user=eamcdani@mail.ubc.ca   
#SBATCH --mail-type ALL  

# paths 
project_path="/tahoma/emsla51366"
reads_path="${project_path}/polishing/nanopore_lrs_trimmed"
racon_x2_path="${project_path}/polishing/racon_x2"
racon_x3_path="${project_path}/polishing/racon_x3"

# polish the 2nd round racon assembly 
singularity exec ../51366-apps.sif racon -m 8 -x -6 -g -8 -w 500 -t 16 ${reads_path}/AD_SIP_trimmed_combined.fastq ${racon_x2_path}/assembly_racon_x2_mapping.sam ${racon_x2_path}/assembly_racon_x2.fasta > ${racon_x3_path}/assembly_racon_x3.fasta

# map to the 2nd racon polished assembly
singularity exec ../51366-apps.sif minimap2 -t 16 -x map-ont -a ${racon_x3_path}/assembly_racon_x3.fasta ${reads_path}/AD_SIP_trimmed_combined.fastq > ${racon_x3_path}/assembly_racon_x3_mapping.sam 