#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=2:0:0
#SBATCH --job-name=anvio-commands
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

# paths 
project_path="/project/6049207/AD_metagenome-Elizabeth"
out_path="${project_path}/re_binning/taxonomy"
assembly_path="${project_path}/re_binning/final_np_polished_contigs_modified.fasta"

# load modules
module load python scipy-stack diamond hmmer prodigal
source ~/projects/rrg-ziels/shared_tools/virtual_envs/anvio_7.0/bin/activate

#generate contigs database
anvi-gen-contigs-database -f ${assembly_path} -o ${out_path}/contigs.db -n 'SAOB_ONT_v2' -T 12

anvi-run-hmms -c ${out_path}/contigs.db --num-threads 12

#get taxonomy
anvi-run-scg-taxonomy -c ${out_path}/contigs.db --num-parallel-processes 12 --num-threads 12 --min-percent-identity 80 --all-hits-output-file scg-taxonomy.txt

anvi-estimate-scg-taxonomy -c ${out_path}/contigs.db -T 12 --metagenome-mode

# export gene calls and scg taxonomy results 
anvi-export-gene-calls -c ${out_path}/contigs.db --gene-caller prodigal -o ${out_path}/gene_calls.txt 

anvi-get-sequences-for-gene-calls -c contigs.db --get-aa-sequences -o saob-ont-proteins.faa

anvi-script-get-hmm-hits-per-gene-call -c contigs.db -o hmm-hits-bacteria.txt --hmm-source Bacteria_71

anvi-script-get-hmm-hits-per-gene-call -c contigs.db -o hmm-hits-archaea.txt --hmm-source Archaea_76



