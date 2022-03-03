#!/bin/bash
#SBATCH --account emsla51366                  
#SBATCH --time 6:0:0 
#SBATCH --nodes 2                      
#SBATCH --ntasks-per-node 24  
#SBATCH --job-name medaka_x1_polishing
#SBATCH --error medaka_x1_polishing.err    
#SBATCH --output medaka_x1_polishing.out          
#SBATCH --mail-user=eamcdani@mail.ubc.ca   
#SBATCH --mail-type ALL  

# paths 
project_path="/tahoma/emsla51366"
reads_path="${project_path}/polishing/nanopore_lrs_trimmed"
racon_x3_path="${project_path}/polishing/racon_x3"
medaka_x1_path="${project_path}/polishing/medaka_x1"
model_path="/tahoma/emsla51366/dbs/r104_e81_sup_g5015_model.tar.gz"

# polish the 3rd round racon assembly with medaka 
singularity exec ${project_path}/51366-apps-new.sif medaka_consensus -i ${reads_path}/AD_SIP_trimmed_combined.fastq -d ${racon_x3_path}/assembly_racon_x3.fasta -o $medaka_x1_path -m ${model_path} -t 24

