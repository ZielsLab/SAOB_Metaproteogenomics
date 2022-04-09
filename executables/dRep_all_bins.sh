#!/bin/bash
#SBATCH --account=def-shallam
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=15G
#SBATCH --time=12:0:0
#SBATCH --job-name=dRep-sakinaw.sh
#SBATCH --output=dRep-sakinaw.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
bins_path="${project_path}/re_binning_all_TS/all_bins"
out_path="${project_path}/re_binning_all_TS/dRep"
dRep_env="/home/eamcdani/virtual_envs/dRep/bin/activate"
bin_stats="${project_path}/re_binning_all_TS/checkM/dRep-input.csv"

# load dRep virtual env
source ${dRep_env}

# load necessary modules 
module load StdEnv/2020 gcc/9.3.0 mash/2.3 mummer fastani prodigal centrifuge
PYTHONPATH=''

# dRep command 
dRep dereplicate ${out_path} -g ${bins_path}/*.fa --genomeInfo ${bin_stats}

deactivate