#!/bin/bash
#SBATCH --account emsla51366
#SBATCH --array=1-13
#SBATCH --time 24:0:0
#SBATCH --nodes 2
#SBATCH --ntasks-per-node 32
#SBATCH --job-name spades_assembly
#SBATCH --error spades_assembly.err
#SBATCH --output spades_assembly.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type ALL

# samples from array 
sample_name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" illumina_samples.txt)

# paths
project_path="/tahoma/emsla51366"
reads_path="${project_path}/illumina_trimmed_filtered"
r1_file="${reads_path}/trim${sample_name}_R1_.fastq"
r2_file="${reads_path}/trim${sample_name}_R2_.fastq"
out_path="${project_path}/assemblies/${sample_name}"

# spades call

singularity exec ${project_path}/51366-apps-new.sif spades.py --meta -t 32 -m 500 -k 21,33,55,77,99,127 -1 $r1_file -2 $r2_file -o ${out_path}