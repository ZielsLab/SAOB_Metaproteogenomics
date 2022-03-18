#!/bin/bash
#SBATCH --account emsla51366
#SBATCH --array=1-15
#SBATCH --time 24:0:0
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --job-name mapMetasToNanopore
#SBATCH --error mapMetasToNanopore.err
#SBATCH --output mapMetasToNanopore.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type ALL

# mapping combo info from array 
mapping_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metagenomeList.txt)
mapping_name=$(basename $mapping_file _R1_.fastq)

# paths
project_path="/tahoma/emsla51366"
reads_path="${project_path}/illumina_trimmed_filtered"
r1_file=$(basename $mapping_file)
r2_file=$(basename $mapping_file _R1_.fastq)_R2_.fastq
r1_path="${reads_path}/${r1_file}"
r2_path="${reads_path}/${r2_file}"
index_path="${project_path}/assemblies/polished_nanopore/bt2/polished_nanopore.fasta"
out_path="${project_path}/mappingResults"
out_name="${mapping_name}-vs-polished_nanopore"

# bowtie2 mapping command 
singularity exec ${project_path}/51366-apps-new.sif bowtie2 -x $index_path -1 $r1_path -2 $r2_path -q --very-sensitive -p 16 -S ${out_path}/${out_name}.sam

# samtools SAM to BAM 
singularity exec ${project_path}/51366-apps-new.sif samtools view -@ 12 -S -b ${out_path}/${out_name}.sam > ${out_path}/${out_name}.bam

# samtools sort and index
singularity exec ${project_path}/51366-apps-new.sif samtools sort -@ 12 ${out_path}/${out_name}.bam -o ${out_path}/${out_name}.sorted.bam
singularity exec ${project_path}/51366-apps-new.sif samtools index ${out_path}/${out_name}.sorted.bam