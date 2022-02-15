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
inpath="${project_path}/05_racon_polish"
outpath="${project_path}/06_medaka_polish"
#model_path="/project/6049207/shared_tools/medaka_models/r104_e81_sup_g5015_model.tar.gz"
model_path="r104_e81_sup_g5015"

#prepare environment
module load gcc python bcftools minimap2 samtools
umask g+w
source ${medaka_path}

#run medaka
medaka_consensus -i ${read_path} \
    -d ${inpath}/assembly_racon_x1.fasta \
    -o $outpath \
    -m ${model_path} \
    -t 24

deactivate
