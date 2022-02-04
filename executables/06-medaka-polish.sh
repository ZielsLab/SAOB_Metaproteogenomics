#!/bin/bash
#SBATCH --account=def-ziels
#SBATCH --gres=gpu:p100l:4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=0
#SBATCH --time=1-12:0
#SBATCH --job-name=06-medaka_polish
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
medaka_path="/project/6049207/shared_tools/virtual_envs/medaka_v1.4.3/bin/activate"
read_path="${project_path}/illumina_raw/MV_MicrOnline_Nov_Dec2021-porechop_trimmed.fastq"
inpath="${project_path}/04_racon_polish"
outpath="${project_path}/05_medaka_polish"

#prepare environment
module load gcc python bcftools minimap2 samtools
umask g+w
source ${medaka_path}

#run medaka
medaka_consensus -i ${read_path} \
	-d ${inpath}/assembly_racon_x1.fasta \
	-o $outpath \
	-m r941_min_sup_g507 \
	-t 24

deactivate