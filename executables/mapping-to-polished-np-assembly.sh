#!/bin/bash
#SBATCH --account=rrg-ziels
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=64G
#SBATCH --time=22:0:0
#SBATCH --job-name=mapping-to-nanopore
#SBATCH --output=%x.out
#SBATCH --mail-user=eamcdani@mail.ubc.ca
#SBATCH --mail-type=ALL

# mapping to the nanopore assembly all illumina short reads and the long nanopore reads to prepare for manual binning in mmgenome and with Anvi'o

#paths
project_path="/project/6049207/AD_metagenome-Elizabeth"
ilm_read_path="${project_path}/illumina_qced"
np_read_path="${project_path}/03_porechop_trimmed/AD_SIP_trimmed_combined.fastq"
assembly_path="${project_path}/re_binning/np_binning_v2_poly/assembly/raconx3_medakax3_raconilmx1_polypolish_contigs.fasta"
out_path="${project_path}/re_binning/np_binning_v2_poly"

# modules
module load StdEnv/2020 gcc/9.3.0 bowtie2 samtools minimap2 metabat/2.14

# map long nanopore reads to nanopore assembly 
minimap2 -t 16 -x map-ont -a ${assembly_path} ${np_read_path} > ${out_path}/mappingResults/final_np_assembly_vs_np_mapping.sam

# create bowtie index of nanopore assembly 
bowtie2-build ${assembly_path} ${out_path}/mappingResults/bt2/np_assembly.fasta 

# map the Sept experimental sample to the nanpore assembly 
bowtie2 -x ${out_path}/mappingResults/bt2/np_assembly.fasta -1 ${project_path}/00_illumina_for_polishing/illumina_qced/R2Sept2020_R1.qced.fastq -2 ${project_path}/00_illumina_for_polishing/illumina_qced/R2Sept2020_R2.qced.fastq -q --very-sensitive -p 16 -S ${out_path}/mappingResults/R2Sept2020.sam

# map all R2 illumina reads to the nanopore assembly
for file in ${ilm_read_path}/R2*_R1.qced.fastq; do
	name=$(basename $file _R1.qced.fastq);
	bowtie2 -x ${out_path}/mappingResults/bt2/np_assembly.fasta -1 ${ilm_read_path}/${name}_R1.qced.fastq -2 ${ilm_read_path}/${name}_R2.qced.fastq -q --very-sensitive -p 16 -S ${out_path}/mappingResults/${name}.sam; 
done 
	
# sort and index files  
for file in ${out_path}/mappingResults/*.sam; do
	name=$(basename $file .sam);
	samtools view -@ 12 -S -b $file > ${out_path}/mappingResults/${name}.bam;
	samtools sort -@ 12 ${out_path}/mappingResults/${name}.bam -o ${out_path}/mappingResults/${name}.sorted.bam;
	samtools index ${out_path}/mappingResults/${name}.sorted.bam;
done 

# get coverage for each sample with jgi_summarize_contigs
jgi_summarize_bam_contig_depths --outputDepth ${out_path}/r2-np-polished-assembly-samples-depth.txt ${out_path}/mappingResults/*.sorted.bam
