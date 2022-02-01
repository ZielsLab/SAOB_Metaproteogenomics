#! /bin/bash

inpath="/var/lib/minknow/data/AD_bioreactor/R2/20220120_1812_MN38981_FAR76126_70891827"
outpath="/home/xlin28/Desktop/AD_Bioreactor/01-basecall_simplex"
guppy_path="/home/xlin28/tools/ont-guppy_6.0.1/bin/guppy_basecaller"

mkdir ${outpath}

$guppy_path -i $inpath -s $outpath \
	-c dna_r10.4_e8.1_sup.cfg \
	-x auto \
	--chunks_per_runner 128 \
	-r