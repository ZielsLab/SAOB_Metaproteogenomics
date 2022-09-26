

#!/usr/bin/env python3

import argparse, os
import warnings 
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np

# Arguments 

parser = argparse.ArgumentParser(description = "Parse TSV from Metapathways to TSV for input into Anvio")
parser.add_argument('mpw_file', metavar='MPW', help='Annotation file from Metapathways in tsv format')
parser.add_argument('--anvio', default='anvio-gene-calls.txt', help="Output: Functional annotation table for input into anvio (Default anvio-gene-calls.txt")
parser.add_argument('--annotation', default='annotation-table.txt', help="Output: Functional annotation table for matching up new gene_callers_id with original locus tag to make anvi'o happy (Default annotation-table.txt")


args = parser.parse_args()

# Input and output files
MPW = args.mpw_file
OUT_ANVIO = open(args.anvio, "w")
OUT_ANNO = open(args.annotation, "w")

# Headers
# OUT_ANNO.write("gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\t\source\tversion\n")

# Parse the TSV file as PD DF
df = pd.read_csv(MPW, sep="\t", names=['contig', 'source', 'type', 'start', 'stop', 'call_type', 'strand', 'partial', 'locus_tag'])
subset_df = df[['locus_tag', 'contig', 'start', 'stop', 'strand', 'partial', 'call_type' 'source']]

new_df = df[['locus_tag', 'contig']]
new_df['start'] = subset_df['start'] - 1
new_df['stop'] = subset_df['stop']
new_df['direction'] = subset_df['strand'].replace(['+', '-'], ['f', 'r'])
new_df['partial'] = subset_df['partial']
new_df['locus_tag'] = new_df['locus_tag'].str.split(r'\s*ID=\s*|\s*\;\s*').str[1]
new_df['source'] = subset_df['source'].str.split(r'\s*:\s*').str[0]
new_df['version'] = subset_df['source'].str.split(r'\s*:\s*').str[1]
new_df['gene_callers_id'] = np.arange(len(new_df))
new_df['source'] = new_df['source'].replace(["Prodigal"], ["prodigal"])
new_df['call_type'] = subset_df['call_type']

final_df = new_df.loc[new_df['source'] == 'prodigal']

anvio_df = final_df[['gene_callers_id', 'contig', 'start', 'stop', 'direction', 'partial', 'call_type', 'source', 'version']]
annotation_df = final_df[['gene_callers_id', 'locus_tag', 'contig', 'start', 'stop', 'direction', 'partial', 'source', 'version']]

anvio_df.to_csv(OUT_ANVIO, sep="\t", index=False)
annotation_df.to_csv(OUT_ANNO, sep="\t", index = False)
