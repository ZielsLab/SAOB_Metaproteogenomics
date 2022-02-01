#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=5:0:0
#SBATCH --job-name=05-racon_polish
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL


#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
read_path="${project_path}/03_porechop_trimmed"
assembly_path="${project_path}/04_flye_assembly_v2"
out_path="${project_path}/05_racon_polish"

#prepare environment
module load singularity
mkdir $out_path

#run racon
singularity exec -B /home -B /project -B /scratch -B /localscratch docker://staphb/racon \
racon -m 8 -x -6 -g -8 -w 500 -t 16 \
$read_path/AD_SIP_trimmed_combined.fastq \
$assembly_path/flye_assembly.sam \
$assembly_path/assembly.fasta > $out_path/assembly_racon_x1.fasta