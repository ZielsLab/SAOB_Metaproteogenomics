#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=25G
#SBATCH --time=5:0:0
#SBATCH --job-name=finalBins_gtdbtk
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL


#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
bins_path="${project_path}/re_binning/combined_bin_set/bins"
out_path="${project_path}/re_binning/combined_bin_set/GTDB"
gtdbtk_path="/home/eamcdani/virtual_envs/gtdbtk/bin/activate"

# modules and gtdbtk env load
source ${gtdbtk_path}
module load StdEnv/2020 gcc/9.3.0 prodigal hmmer pplacer fastani fasttree mash/2.3

# export gtdbtk data path
PYTHONPATH=''
export GTDBTK_DATA_PATH=/home/eamcdani/projects/rrg-ziels/shared_tools/databases/gtdb/release202

# gtdbtk command 

gtdbtk classify_wf --cpus 16 --extension fa --genome_dir ${bins_path}/ --out_dir ${out_path}