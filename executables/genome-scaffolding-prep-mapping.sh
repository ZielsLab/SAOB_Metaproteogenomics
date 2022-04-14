#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12G
#SBATCH --time=24:0:0
#SBATCH --job-name=bin8.2-archaea-mapping
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
read_path="${project_path}/00_illumina_for_polishing/illumina_qced"
ont_read_path="${project_path}/03_porechop_trimmed"
assembly_path="${project_path}/re_binning/np_binning/scaffolding"
out_path="${project_path}/re_binning/np_binning/scaffolding"

# load modules
module load bowtie2 minimap2/2.17

# build bowtie2 index
bowtie2-build ${assembly_path}/bin8.2.fa ${out_path}/bin8.2.fa

# map Illumina reads with bowtie2 
bowtie2 -x ${out_path}/bin8.2.fa -1 ${read_path}/R2Sept2020_R1.qced.fastq -2 ${read_path}/R2Sept2020_R1.qced.fastq -q --no-unal --very-sensitive -p 8 -S ${out_path}/bin8.2.fa_ilm.sam

# map Nanopore reads with minimap2 
minimap2 -t 8 -x map-ont -a ${assembly_path}/bin8.2.fa ${ont_read_path}/AD_SIP_trimmed_combined.fastq > ${out_path}/bin8.2.fa_np.sam 
