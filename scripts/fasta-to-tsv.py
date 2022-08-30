import argparse, os
from Bio import SeqIO 

parser = argparse.ArgumentParser(description = "Convert fasta to tsv")
parser.add_argument("fasta_file", metavar="FASTA", help="Fasta file input")
parser.add_argument("output", default="proteins-table.tsv", help="Output: tsv table of multi-fasta file")

args = parser.parse_args()

FASTA = args.fasta_file 
OUT_TSV = open(args.output, "w") 

# header 
OUT_TSV.write("locus_tag\tAA_sequence\n")

for record in SeqIO.parse(FASTA, "fasta"): 
    locus_tag, aa = (record.description, record.seq)
    OUT_TSV.write('%s\t%s\n' % (locus_tag, aa))