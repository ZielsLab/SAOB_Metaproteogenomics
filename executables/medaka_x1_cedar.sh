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
read_path="${project_path}/03_porechop_trimmed/AD_SIP_trimmed_combined.fastq"
medaka_path="/project/6049207/shared_tools/virtual_envs/medaka_v1.5/bin/activate"
inpath="${project_path}/re_polishing/racon_x3"
outpath="${project_path}/re_polishing/medaka_x1"
#model_path="/project/6049207/shared_tools/medaka_models/r104_e81_sup_g5015_model.tar.gz"
model_path="/home/eamcdani/projects/rrg-ziels/shared_tools/medaka_models/r104_e81_sup_g5015_model.tar.gz"

#prepare environment
module load gcc python bcftools minimap2 samtools
umask g+w
source ${medaka_path}

#run medaka
medaka_consensus -i ${read_path} \
	-d ${inpath}/assembly_racon_x3.fasta \
	-o $outpath \
	-m ${model_path} \
	-t 24

deactivate