#!/bin/bash
#SBATCH --account emsla51366                  
#SBATCH --time 24:0:0 
#SBATCH --nodes 2                      
#SBATCH --ntasks-per-node 32  
#SBATCH --job-name spades_assembly
#SBATCH --error spades_assembly.err    
#SBATCH --output spades_assembly.out          
#SBATCH --mail-user=eamcdani@mail.ubc.ca   
#SBATCH --mail-type ALL  

# paths
project_path="/tahoma/emsla51366"
reads_path="${project_path}/illumina_trimmed_filtered"
out_path="${project_path}/assemblies/R2Sept2020"

# spades call 

singularity exec ${project_path}/51366-apps-new.sif spades.py --meta -t 32 -m 500 -k 21,33,55,77,99,127 -1 ${reads_path}/trimR2Sept2020_R1_.fastq -2 ${reads_path}/trimR2Sept2020_R2_.fastq -o ${out_path}
