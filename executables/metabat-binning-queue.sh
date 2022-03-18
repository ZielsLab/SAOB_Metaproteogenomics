#!/bin/bash
#SBATCH --account emsla51366
#SBATCH --array=1-15
#SBATCH --time 12:0:0
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --job-name SAOB_binning
#SBATCH --error SAOB_binning.err
#SBATCH --output SAOB_binning.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type ALL

# assemblies from array
assembly_name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" saob_assemblies.txt)

# paths
project_path="/tahoma/emsla51366"
mapping_dir="${project_path}/mappingResults"
out_path="${project_path}/binningResults"
assembly_path="${project_path}/assemblies/${assembly_name}/contigs.fasta"

# binning commands 
singularity exec ${project_path}/51366-apps-new.sif jgi_summarize_bam_contig_depths --outputDepth ${out_path}/${assembly_name}-depth.txt ${mapping_dir}/${assembly_name}*.sorted.bam

singularity exec ${project_path}/51366-apps-new.sif metabat2 -i $assembly_path -a ${out_path}/${assembly_name}-depth.txt -o ${out_path}/${assembly_name}/${assembly_name}-bin
