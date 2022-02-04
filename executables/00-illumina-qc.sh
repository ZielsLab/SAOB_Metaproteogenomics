#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=20G
#SBATCH --time=1-12:0
#SBATCH --job-name=00-illumina-qc
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
read_path="${project_path}/illumina_raw"
out_path="${project_path}/illumina_qced"

# prepare environment
module load bbmap
mkdir $out_path

# bbduk QC in bbmap suite
bbduk.sh in=${read_path}/R2Sept2020_R1_.fastq in2=${read_path}/R2Sept2020_R2_.fastq out=${out_path}/R2Sept2020_R1.qced.fastq out2=${out_path}/R2Sept2020_R2.qced.fastq qtrim=r trimq=10 maq=10

