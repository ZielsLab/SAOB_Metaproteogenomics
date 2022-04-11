#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=64G
#SBATCH --time=22:0:0
#SBATCH --job-name=rebinning_checkM_manual
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
bins_path="${project_path}/re_binning/np_binning/manual_binning/fastas"
out_path="${project_path}/re_binning/np_binning/manual_binning/checkM"
checkm_path="/home/eamcdani/projects/rrg-ziels/shared_tools/virtual_envs/checkM/bin/activate"

#prepare environment
module load gcc python hmmer prodigal pplacer
source ${checkm_path}

# checkM commands
checkm lineage_wf -x .fa -t 8 ${bins_path}/ ${out_path}
cd ${out_path}
checkm qa lineage.ms -o 2 -f checkm.out ./ --tab_table

# format out file for table and dRep input
awk -F "\t" '{print $1"\t"$2"\t"$6"\t"$7"\t"$8"\t"$9"\t"$12"\t"$19}' checkm.out > checkm_stats.tsv
awk -F "\t" '{print $1".fa,"$3","$4}' checkm_stats.tsv > dRep-input.csv