

#!/usr/bin/env python

import argparse, os
import warnings 
warnings.filterwarnings("ignore")
import pandas as pd

# Arguments 

parser = argparse.ArgumentParser(description = "Parse TSV from Metapathways to TSV for input into Anvio")
parser.add_argument('mpw_file', metavar='MPW', help='Annotation file from Metapathways in tsv format')
parser.add_argument('--annotation', default='anvio_gene_calls.txt', help="Output: Functional annotation table for input into anvio (Default anvio_gene_calls.txt")

args = parser.parse_args()

# Input and output files
MPW = args.mpw_file
OUT_ANNO = open(args.annotation, "w")

# Headers
# OUT_ANNO.write("gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tcall_type\t\source\tversion\n")

# Parse the TSV file as PD DF
df = pd.read_csv(MPW, sep="\t", names=['contig', 'source', 'type', 'start', 'stop', 'call_type', 'strand', 'partial', 'locus_tag'])
subset_df = df[['locus_tag', 'contig', 'start', 'stop', 'strand', 'partial', 'call_type', 'source']]

new_df = df[['locus_tag', 'contig']]
new_df['start'] = subset_df['start'] - 1
new_df['stop'] = subset_df['stop']
new_df['direction'] = subset_df['strand'].replace(['+', '-'], ['f', 'r'])
new_df['partial'] = subset_df['partial']
new_df['call_type'] = subset_df['call_type'].replace(['.'], ['1'])
new_df['gene_callers_id'] = new_df['locus_tag'].str.split(r'\s*ID=\s*|\s*\;\s*').str[1]
new_df['source'] = subset_df['source'].str.split(r'\s*:\s*').str[0]
new_df['version'] = subset_df['source'].str.split(r'\s*:\s*').str[1]

final_df = new_df[['gene_callers_id', 'contig', 'start', 'stop', 'direction', 'partial', 'call_type', 'source', 'version']]

final_df.to_csv(OUT_ANNO, sep="\t", index = False)
