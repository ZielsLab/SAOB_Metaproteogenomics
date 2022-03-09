#!/bin/bash
#SBATCH --account emsla51366
#SBATCH --array=1-14
#SBATCH --time 6:0:0
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --job-name bowtie_index_assemblies
#SBATCH --error bowtie_index_assemblies.err
#SBATCH --output bowtie_index_assemblies.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type ALL

# samples from array 
sample_name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bioreactor-samples.txt)

# paths
project_path="/tahoma/emsla51366"
assembly_path="${project_path}/assemblies/${sample_name}"
out_path="${assembly_path}/bt2"

# make bt2 directory
mkdir ${out_path}

# bowtie2 index
singularity exec ${project_path}/51366-apps-new.sif bowtie2-build ${assembly_path}/contigs.fasta ${out_path}/${sample_name}.fasta