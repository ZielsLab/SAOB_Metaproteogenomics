#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --array=1-15
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=50G
#SBATCH --time=12:0:0
#SBATCH --job-name=inStrain-quick-profile-saob
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

# array jobs 
mapping_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" saob_sorted_bams.txt)

# sample name
sample_name=$(basename $mapping_file .sorted.bam)

#paths
inStrain_env="/home/eamcdani/virtual_envs/inStrain/bin/activate"
project_path="/project/6049207/AD_metagenome-Elizabeth"
bam_path="${project_path}/re_binning/combined_bin_set/mappingResults"
fasta="${project_path}/re_binning/combined_bin_set/mappingResults/bins/all_SAOB_bins.fasta"
out_path="${project_path}/re_binning/combined_bin_set/inStrain/quick_profiles"
out_name="${out_path}/${sample_name}"
stb="${project_path}/re_binning/combined_bin_set/mappingResults/bins/scaffolds-to-bins.tsv"


# load modules and programs
source ${inStrain_env}
export PATH="/home/eamcdani/.cargo/bin/:$PATH"
module load samtools

# cd to mapping directory 
cd ${bam_path}

# quick profile command
inStrain quick_profile -p 2 -s $stb -o $out_name $mapping_file $fasta

# cd to the results directory to change the genomeCoverage filename
cd ${out_name}
mv genomeCoverage.csv ${sample_name}-genomeCoverage.csv
