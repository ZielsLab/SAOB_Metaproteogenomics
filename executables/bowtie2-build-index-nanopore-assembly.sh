#!/bin/bash
#SBATCH --account emsla51366                  
#SBATCH --time 5:0:0 
#SBATCH --nodes 1                      
#SBATCH --ntasks-per-node 16  
#SBATCH --job-name build_np_assemb_index
#SBATCH --error build_np_assemb_index.err    
#SBATCH --output build_np_assemb_index.out          
#SBATCH --mail-user=eamcdani@mail.ubc.ca   
#SBATCH --mail-type ALL  

# paths
project_path="/tahoma/emsla51366"
assembly_path="${project_path}/assemblies/polished_nanopore"
out_path="${assembly_path}/bt2"

# make bt2 directory
mkdir ${out_path}

# bowtie2 index
singularity exec ${project_path}/51366-apps-new.sif bowtie2-build ${assembly_path}/contigs_medaka.consensus_racon.ilm_x1.fasta ${out_path}/polished_nanopore.fasta
