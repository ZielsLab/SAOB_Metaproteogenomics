# Methanothermobacter Pangenomics with Anvio 

This page documents steps for creating a pangeome of Methanothermobacter references and genomes assembled in this study using Anvi'o. The ultimate goal of this analysis is to identify core and accessory sets of groups of proteins to then look at the protein expression dyanmics of these groups of proteins of the Methanothermobacter in the SIP microcosm experiment. This workflow follows the Anvi'o pangenomics workflow: https://merenlab.org/2016/11/08/pangenomics-v2/ 

## Collecting reference genomes

First I searched the GTDB for "Methanothermobacter" and downloaded those genomes with `ncbi-genome-download`. Then I used only genomes with above 90% completeness and less than 10% redundancy, for a total of 24 reference genomes, which are listed in the `methanothermobacter_gtdb_metadata.csv` file. From this study I assembled 4 different Methanothermobacter genomes, two of which `METHANO1` and `METHANO2` are highly enriched in lab-scale bioreactors along with a putative syntrophic acetate oxidizing bacteria `DTU030_1`. I only included three of the four Methanothermobacter genomes in this analysis because the _Methanothermobacter thermoautotrophicus_ genome only has 79% completion and is 229 contigs, and therefore isn't that great of a genome for this analysis. Whereas the other 3 Methanothermobacter genomes were assembled with Nanopore data, have decently high completion, and are on 5 or 6 contigs - so pretty good genomes. 

## Annotation of Methanothermobacter references 

For the 3 Methanothermobacter genomes assembled in this study, I used the MetaPathways2 program for functional annotation and exploring pathways - which does this by calling genes with Prodigal and making functional annotations with Prokka. Since I have to keep these locus tags the same to compare to in the proteomics data, I am keeping these annotations and will have to load these into Anvi'o with `external-gene-calls`. However for the other reference genomes, I can load these into Anvi'o and have the pipeline itself perform annotation. 

## Prepping Genomes for Anvi'o

We will process the reference genomes and 3 Methanothermobacter SAOB study genomes separately for the first few steps since the annotation for the 3 study genomes is different. 

First reformat the contig names of the reference genomes with `anvi-script-reformat-fastas`: 

```
for file in *.fna; 
    do name=$(basename $file .fna); 
    anvi-script-reformat-fasta $file -o $name-reformatted.fasta --simplify-names; 
done
```

Then create contigs databases for each reformatted fasta file: 

```
for file in *-reformatted.fasta; do
    name=$(basename $file .1-reformatted.fasta);
    anvi-gen-contigs-database -f $file -o ../contigs_dbs/$name.db;
done
```

The GFF files from Metapathways from the 3 Methanothermobacter SAOB study genomes are then parsed with the `metapathways-tsv-to-anvio.py` script with:

```
for file in *.annot.gff; 
    do name=$(basename $file .annot.gff); 
    python3 ../../../scripts/metapathways-tsv-to-anvio.py $file --annotation $name-anvio-table.tsv; 
done
```

These can then be created into contigs databases with the `--external-gene-calls` flag: 

```
for file in *.fasta; 
    do name=$(basename $file .fasta);
    anvi-gen-contigs-database -f $file -o ../contigs_dbs/$name.db --external-gene-calls $name-anvio-table.tsv;
done
```

### Functional Annotation 

To add functional annotation information, we do this to each contigs database. The different options you can populate a contigs database with are [described here](https://anvio.org/help/main/artifacts/contigs-db/). For this pangenomics analysis, we will populate each contigs database with the HMM collections for archaea and the NCBI COGs annotations.

First run the archaeal HMMs (doesn't make sense to run all the others because these are all Methanothermobacter genomes)

```
for file in *.db; do
    anvi-run-hmms -c $file -I Archaea_76; 
done
```

Then add the NCBI COGs annotations: 

```
for file in *.db; do 
    anvi-run-ncbi-cogs -c $file -T 6; 
done
```

Then create an anvi'o genomes storage for the entire collection of genomes. Make a list of all the contigs databases with `for file in *.db; do name=$(basename $file .db); echo -e $name"\t"$file; done > methano_genomes_storage.txt`, append with colum names `name` and `contigs_db_path`, and then create the genomes storage: `anvi-gen-genomes-storage -e methano_genomes_storage.txt -o METHANO_GENOMES.db`


## Running the Pangenome Analysis

```
anvi-pan-genome -g METHANO_GENOMES.db \
                --project-name "Methanothermobacter_Pan" \
                --output-dir METHANO \
                --num-threads 6 \
                --minbit 0.5 \
                --mcl-inflation 2
```

This will use DIAMOND to perform all-v-all comparisons and an MCL inflation score of 2. The higher the MCL inflation score (such as 10) is for very closely related strains, and the lower MCL inflation score of 2 is for more distantly related groups of genomes. 

Then compute average nucleotide identity (ANI) between all genomes:

```
anvi-compute-genome-similarity --external-genomes methano_genomes_storage.txt \
                               --program pyANI \
                               --output-dir ANI \
                               --num-threads 6 \
                               --pan-db METHANO/Methanothermobacter_Pan-PAN.db
```    

This will add a layer to the pangenomics visualization output to view pairwise ANI comparisons between all genomes, instead of having to do this separatley outside of Anvi'o. 

## Interactively Visualize the Pangenome 

If you are on a remote compute cluster, [this tutorial](https://merenlab.org/2015/11/28/visualizing-from-a-server/) explains how to run an SSH tunnel to view the interactive interface. 

Once you are ssh'ed into the cluster, view the pangenome with: 

```
anvi-display-pan -p METHANO/Methanothermobacter_Pan-PAN.db -g METHANO_GENOMES.db --server-only -P 8080
```

Then start your browser on your local computer and type the address `http://localhost:8080` 

### 