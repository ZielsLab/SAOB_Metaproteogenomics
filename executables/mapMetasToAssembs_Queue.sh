#!/bin/bash
#SBATCH --account emsla51366
#SBATCH --array=1-210
#SBATCH --time 24:0:0
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --job-name mapMetasToAssembsQueue
#SBATCH --error mapMetasToAssembsQueue.err
#SBATCH --output mapMetasToAssembsQueue.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type ALL

# mapping combo info from array 
sample_name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" illumina_sample_mappingCombos.txt | awk -F "\t" '{print $1}')

mapping_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" illumina_sample_mappingCombos.txt | awk -F "\t" '{print $2}')

# paths
project_path="/tahoma/emsla51366"
reads_path="${project_path}/illumina_trimmed_filtered"
mapping_name=$(basename $mapping_file _R1_.fastq)
r1_file=$(basename $mapping_file)
r2_file=$(basename $mapping_file _R1_.fastq)_R2_.fastq
r1_path="${reads_path}/${r1_file}"
r2_path="${reads_path}/${r2_file}"
index_path="${project_path}/assemblies/${sample_name}/bt2/${sample_name}.fasta"
out_path="${project_path}/mappingResults"
out_name="${mapping_name}-vs-${sample_name}"

# bowtie2 mapping command 
singularity exec ${project_path}/51366-apps-new.sif bowtie2 -x $index_path -1 $r1_path -2 $r2_path -q --very-sensitive -p 16 -S ${out_path}/${out_name}.sam

# samtools SAM to BAM 
singularity exec ${project_path}/51366-apps-new.sif samtools view -@ 12 -S -b ${out_path}/${out_name}.sam > ${out_path}/${out_name}.bam

# samtools sort and index
singularity exec ${project_path}/51366-apps-new.sif samtools sort -@ 12 ${out_path}/${out_name}.bam -o ${out_path}/${out_name}.sorted.bam
singularity exec ${project_path}/51366-apps-new.sif samtools index ${out_path}/${out_name}.sorted.bam
