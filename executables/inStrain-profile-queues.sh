#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --array=1-343
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=12:0:0
#SBATCH --job-name=inStrain-profiles-queue
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

# array jobs 
bam_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" inStrain-profiles-queue.txt | awk -F "\t" '{print $1}')
genome=$(sed -n "${SLURM_ARRAY_TASK_ID}p" inStrain-profiles-queue.txt | awk -F "\t" '{print $2}')

# sample name
sample_name=$(basename $bam_file .sorted.bam)
genome_name=$(basename $genome .fa)

#paths
inStrain_env="/home/eamcdani/virtual_envs/inStrain/bin/activate"
project_path="/project/6049207/AD_metagenome-Elizabeth"
bam_path="${project_path}/re_binning/combined_bin_set/mappingResults"
mapping="${bam_path}/${bam_file}"
fasta_path="${project_path}/re_binning/combined_bin_set/mappingResults/bins/"
fasta="${fasta_path}/${genome}"
genes="${fasta_path}/${genome_name}.genes.fna"
out_path="${project_path}/re_binning/combined_bin_set/inStrain/"
out_name="${out_path}/${sample_name}-vs-${genome_name}"
stb="${project_path}/re_binning/combined_bin_set/mappingResults/bins/scaffolds-to-bins.tsv"

# load modules and programs
source ${inStrain_env}
export PATH="/home/eamcdani/.cargo/bin/:$PATH"
module load samtools

# cd to mapping directory 
cd ${bam_path}

# quick profile command
inStrain profile $mapping $fasta -o $out_name -p 8 -g $genes -s $stb

