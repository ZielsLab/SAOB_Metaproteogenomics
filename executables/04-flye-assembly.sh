#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=150G
#SBATCH --time=1-12:0
#SBATCH --job-name=04-flye-assembly
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL


#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
flye_path="/project/6049207/shared_tools/conda_envs/flye"
conda_path="/home/eamcdani/miniconda3/bin"
inpath="${project_path}/03_porechop_trimmed"
outpath="${project_path}/04_flye_assembly"

#prepare environment
source ${conda_path}/activate ${flye_path}
mkdir ${outpath}

#run flye
flye --nano-hq ${inpath}/AD_SIP_trimmed_combined.fastq -o $outpath -t 16 --meta --scaffold

#get mapping stats
module load minimap2/2.17

minimap2 -t 16 -x map-ont -a ${outpath}/assembly.fasta \
        ${inpath}/AD_SIP_trimmed_combined.fastq > ${outpath}/flye_assembly.sam

source deactivate