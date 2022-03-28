#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=12G
#SBATCH --time=14:0:0
#SBATCH --job-name=taxonomy-prep
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
out_path="${project_path}/re_binning/taxonomy"
assembly_path="${project_path}/re_binning/final_np_polished_contigs_modified.fasta"
db="/project/6049207/shared_tools/databases/diamond/annotree.dmnd"

module load diamond

diamond blastx -d ${db} -q ${assembly_path} -o ${out_path}/contigs-annotree-tax.daa -f 100 --long-reads -p 12

#transfer that daa file to EAM computer to run megan6
# /Applications/MEGAN/tools/daa-meganizer -i contigs-annotree-tax.daa -mdb megan-mapping-annotree-June-2021.db -lg
#  /Applications/MEGAN/tools/daa2info -i contigs-annotree-tax.daa -o contigs_megan_tax.tsv -l -m -r2c GTDB -p true -r true