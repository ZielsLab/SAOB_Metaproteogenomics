#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=24:0:0
#SBATCH --job-name=kofam-annotate-methanothermobacter.sh
#SBATCH --output=kofam-annotate-methanothermobacter.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

# run all genome proteins concatenated together with kofamkoala - runs faster than protein FASTAs split individually for each genome
# paths 
#paths
fasta_path="/home/eamcdani/scratch"
kofam_path="/home/eamcdani/bin/kofam_scan-1.3.0"

# load modules 
module load ruby/2.7.1 hmmer 

${kofam_path}/exec_annotation ${fasta_path}/all-methanothermobacter-proteins.faa -o /home/eamcdani/scratch/all-methanothermobacter-kofam-annotations.txt -p ${kofam_path}/profiles/ -k ${kofam_path}/ko_list --cpu 8

# after annotate, modify the output
grep -w '*' /home/eamcdani/scracth/all-methanothermobacter-kofam-annotations.txt | awk -F " " '{print $2"\t"$3}' | sed 's/_/-/g' | sed 's/~/_/g' > /home/eamcdani/scratch/all-methanothermobacter-kofam-modified.txt