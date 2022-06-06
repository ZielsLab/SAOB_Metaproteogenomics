#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --array=1-3
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=60G
#SBATCH --time=48:0:0
#SBATCH --job-name=unicycler-hybrid-assembly-scratch
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
dirs_path="${project_path}/re_binning/np_binning_v2_poly/unicycler_hybrid"
scratch_path="/home/eamcdani/scratch"

# load modules
module load singularity 

# assembly directory name based from the array job
dir_name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" unicycler_dir_files.txt | awk -F "\t" '{print $1}')

# cd into the directory 
cd ${dirs_path}/${dir_name}

# unicycler command 
singularity exec -B /home -B /project -B /scratch -B /localscratch ${project_path}/singularity_images/unicycler_v0.4.9.sif unicycler -1 ilm_assembly_vs_ilm_reads_1.fastq -2 ilm_assembly_vs_ilm_reads_2.fastq -l np_assembly_vs_long_reads.fastq -o ${scratch_path}/${dir_name}_unicycler