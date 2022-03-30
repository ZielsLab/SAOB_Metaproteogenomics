#!/bin/bash
#SBATCH --account emsla51366                  
#SBATCH --time 48:0:0 
#SBATCH --nodes 10                      
#SBATCH --ntasks-per-node 32  
#SBATCH --job-name inoculum_spades_assembly
#SBATCH --error inoculum_spades_assembly.err    
#SBATCH --output inoculum_spades_assembly.out          
#SBATCH --mail-user=eamcdani@mail.ubc.ca   
#SBATCH --mail-type ALL  

# paths
project_path="/tahoma/emsla51366"
reads_path="${project_path}/illumina_trimmed_filtered"
out_path="${project_path}/assemblies/InoculumNov2019"

# spades call 

singularity exec ${project_path}/51366-apps-new.sif spades.py --meta -t 100 -m 2000 -k 21,33,55,77,99,127 -1 ${reads_path}/trimInoculumNov2019_R1_.fastq -2 ${reads_path}/trimInoculumNov2019_R2_.fastq -o ${out_path}
