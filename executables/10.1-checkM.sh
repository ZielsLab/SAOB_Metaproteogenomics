#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=64G
#SBATCH --time=22:0:0
#SBATCH --job-name=10.1-checkM
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
bins_path="${project_path}/10_binning/bins"
out_path="${project_path}/10_binning/checkM"
checkm_path="/home/eamcdani/projects/rrg-ziels/shared_tools/virtual_envs/checkM/bin/activate"

#prepare environment
module load gcc python hmmer prodigal pplacer
source ${checkm_path}

# checkM commands
checkm lineage_wf -x .fa -t 8 ${bins_path}/ ${out_path}
checkm qa ${out_path}/lineage.ms -o 2 -f ${out_path}/checkm.out --tab_table