#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=150GB
#SBATCH --time=1-0:0
#SBATCH --job-name=02-adapter_trimming
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

umask g+w

#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
artic_path="/project/6049207/shared_tools/artic"
conda_path="/home/eamcdani/miniconda3/bin"
inpath="${project_path}/02-split-run_8"
outpath="${project_path}/03_porechop_trimmed"

#prepare environment
source ${conda_path}/activate ${artic_path}

mkdir ${outpath}

for i in {1..10}
do
porechop -i $inpath/split_$i \
	-o $outpath/AD_SIP-trimmed_${i}.fastq \
	--threads 16 --adapter_threshold 95.0
done

deactivate