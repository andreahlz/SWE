# SWE Project Logbook

## NEXT TO DO:
- testing pipline
- long_reads and short_reads nor fair comparison possible (at the moment): 
    - short_reads: number of reads is defined (--n_reads), therefore independent of abbundance and genome size
    - long_reads: number of reads is coverage based --quantity "${coverage}x", therefore dependent of genome size and abundance
    -> idenifying which approach is correct and then using it for both processes

    coverage, not equal read count?
        
    Fixed read count creates unequal information content:
    - 20 short reads (300 bp) = 6,000 bp total information
    - 20 long reads (10,000 bp) = 200,000 bp total information
    - Long reads would have 33× more information (biased comparison)

    Equal coverage ensures equal genomic information:
    - Both: 10× coverage
    - Short reads: ~167,000 reads × 300 bp = 50 Mb
    - Long reads: ~5,000 reads × 10,000 bp = 50 Mb
    - Same information content, different fragmentation

    -> NEW PROBLEM: Different number of embeddings per read_type or is it not relevant? Isnt analysis based on quality not quantity
    
- testing script for simulation check: abundance correct? stats on number of reads and from which species
- which error model should we use?
- nextlow config with parameters for slurm 
- check how to run things on lisc 
-  add labels/species 




### update: pull request changes von 08.01.2026

files changed:
- enviroment.yml
- simulate_short_reads.sh
- .gitignore

new files:
- config.env
- test_abbundance.txt 
- multifasta\test_mock_gut_genomes.fna

info on changes:
- enviroment.yml
    enviroment.yml erstellt auf LISC deshalb pakete linux-spezifisch -> funktioniert nur auf linux enviroment. für windows nur nutzung mit zb wsl
    tasuch von reinfolge von channels wiel conda-forge größer ist zuerst also conda-forge 20000 packages und bioconda 8000 -> reihnfolge nach prioriät und großre und conda-forge hat oft neuere versionen als bioconda
    wenn conda-forge zuerst ist kriegt man dann meist die neusten versionen
    wenn nutzung von wsl ist wahrscheinlich conda nicht installiert auf WSL/Ubuntu subsystem -> deshalb muss installiert werden. keine ahnung aber scheinbar ist in dem fall miniconda besser als anaconda. ich hab das mit den befehlen innerhalb wsl gemacht:
    cd ~
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh

    dann termianl schliesßen und neues öffnene controlle ob instaliiert mit version check:
    conda –version
    git ist wharscheinlich nicht uptodate -> bringen auf selben stand wie github repo

- .gitignore:
    fasta files und fastq files ausgeschlossen aus git repo

- config.env:
    relative paths und parameter hier gespeichert

- simulate_short_reads.sh:
    integration von config file

- test_abbundance.txt und multifasta\test_mock_gut_genomes.fna:
    um zu sehen das simulation functioniert
    simulation funktioniert 
    zwei fastq files werden erstellt



### update: pull request changes von 13.12.2025

files changed:
- config.env

new files:
- simulate_long_reads.sh

info on changes:
- config.env:
  added paragraph: #badread parameters


Added to LiSC manually:
- genome_multifasta="data/genomes/multifasta/test_mock_gut_genomes.fna"
- abundance_file="data/abundances/test_abbundance.txt"
- long_reads_file = "data/scripts/simulate_long_reads.sh"
- "config-env"
Manually because of git permission issues 


### UPDATE: andrea/combine-embeddings

main changes:
- Embedding combination process
- K-means clustering on combined embeddings
- t-SNE visualization of clustering results

files changed:
- scripts/embedding_calculation_paralellized.nf
- dev_doc.md

new files:
- scripts/combined_embeddings.py
- scripts/kmeans_clustering.py
- scripts/visualize_embeddings.py

scripts/combined_embeddings.py: 
Purpose: Combines individual species embedding CSVs into single files
Input: Individual CSV files from calculate_embedding_for_tsv.py (per species)
Output: 
- combined_embeddings.csv - All embeddings combined (for K-means)
- combined_embeddings_stand.csv - Standardized embeddings (for t-SNE)

scripts/kmeans_clustering.py
Purpose: Perform K-means clustering on combined embeddings
Input: combined_embeddings.csv
Output:
- kmeans_results.csv - True labels + Predicted cluster IDs
- kmeans_results_metrics.txt - Clustering quality metrics
script based on eval_clustering_classification_changed.py main diff is adaption of funtions for nextflow pipline
short summary:
- K-means algorithm & parameters
- Multiple random seeds (0-4) for robustness
- Metrics calculation: Purity, Completeness, ARI, NMI
- Label preparation (string → numeric conversion)
- CSV input/output

scripts/visualize_embeddings.py
Purpose: t-SNE visualization of embeddings with clustering results
Input:
combined_embeddings_stand.csv - Standardized embeddings
kmeans_results.csv - Clustering results
Output:
embeddings_visualization.png - Three-panel plot

scripts/embedding_calculation_paralellized.nf
Changes: Extended workflow with three new processes

Current Nextflow Pipline
TSV files
-> 
Calculate embeddings (parallel per species)
-> 
Combine embeddings 
-> 
K-means clustering 
-> 
t-SNE visualization 
-> 
PNG output + CSV results + metrics file

### update: pull request changes von 25.01.2025
short summary: previous testing pipelines did not for set abbundances for each species. New pipeline only needs an abbundance.txt as input and outputs everthing from downloaded genomes, simulated reads, calculated embeddings to visualization automaticaly in one continues pipeline. currently only for short reads.

input format: txt
(base) andyfely@Andys:/mnt/c/Users/andyf/SWE/data/input$ head input.txt
Accession_ID Abundance
GCF_000025985.1 0.4

-> ncbi accesion_ids are used to download Genome sequences.
(script: scripts\download_genome.py)
(base) andyfely@Andys:/mnt/c/Users/andyf/SWE/data/genomes/fasta_files$ head GCF_000025985.1.fna 
>NC_003228.3 Bacteroides fragilis NCTC 9343, complete sequence
TTATCAACACCTATGTTAACAAGAAAAGAATTACTTTTGCAACATACTAACAGAAACGACATCATCATGCGAAAATTGAA
AATAACCGAGCTGAACCGGATAAGTATAGAAGAGTTTAAAGAAGCTGATAAATTGCCTTTAGTTGTAGTGTTGGACGATA
TACGGAGTTTGCATAATATCGGTTCTGTGTTTCGTACGGCAGATGCTTTCCGGATTGAATGTATTTATCTGTGTGGAATT

since some genome fasta have multiple sequences (f.e. chromosome 1 and 2) with multiple header, the faster files are processed so only 1 seqeunce with 1 header per genome fasta and then comined to 1 multifasta. Additionally a csv is created where the informations on all included species are stored.
(script: scripts\)
(base) andyfely@Andys:/mnt/c/Users/andyf/SWE/data/genomes$ head header_info.csv 
accession,headers,abundance
GCF_000025985.1,"NC_003228.3 Bacteroides fragilis NCTC 9343, complete sequence; NC_006873.1 Bacteroides fragilis NCTC 9343 plasmid pBF9343, complete sequence",0.4
GCF_000196555.1,"NC_015067.1 Bifidobacterium longum subsp. longum JCM 1217, complete sequence",0.05
(base) andyfely@Andys:/mnt/c/Users/andyf/SWE/data/genomes$ head mock_test.fna 
>GCF_000025985.1
TTATCAACACCTATGTTAACAAGAAAAGAATTACTTTTGCAACATACTAACAGAAACGACATCATCATGCGAAAATTGAA
AATAACCGAGCTGAACCGGATAAGTATAGAAGAGTTTAAAGAAGCTGATAAATTGCCTTTAGTTGTAGTGTTGGACGATA

read simualtion for short reads:
(script: scripts\)
(base) andyfely@Andys:/mnt/c/Users/andyf/SWE/data/short_reads/simulated_reads$ head mock_test.fastq 
@GCF_000025985.1_0_0/1
TTCTTGATTTGTATCCATGAAAAAAGGGCATAGCGATTATCGCTATGCCCTGCAAATATATTCATTTTATTTATAACGACAATACTGATTGTATATCTACGTCTTGGTCGTTTTCAAATACTTGTTTTCCGTTCAGTTCCTTTAACTCAATGATGAAGTTTACGTATACCTTTTTCGGATGGAGCTCTTTCACGAGGTTGCAGGCGGCTTTCATAGTACCTCCGGTAGCCAGTAAGTCATCGTGAGACAATACGACGTCGTTCTCGTTCAATGCATCTTTGTGTATCTGCACGGTGTCTTT
+
CCCCC@GG@FF9EC@GFCGEGFC@GFFFGFGG@EEA@FCFGGGFEG@F@GFA97GGFGGCGGGFF,GGECFF@GGGCGFF<DCG@6FGGGG?GGFFEGFF:=FGGGF,FFGDAB@GGDFEAGFB<CFDGGG8FF@FEC<E4G:CBGEDDACG<,FFBB=FD:GC9,EF?FFAC<:F,E:?5+=<4+9@<EEF,3@9**,6*FBE7B*3*F3*4G<=7G:DC,E9/+0,43>+6+?1,+9*2*><*+*=<C9)2A4C7++=5)/*1*@+:*++0*;*9*F*8:*1@C:07:D)++2*).0)>

convertion from fastq to tsv for dnaberts input
(script: scripts\)
(base) andyfely@Andys:/mnt/c/Users/andyf/SWE/data/short_reads/simulated_reads$ head mock_test.tsv 
sequence        bin_id
TTCTTGATTTGTATCCATGAAAAAAGGGCATAGCGATTATCGCTATGCCCTGCAAATATATTCATTTTATTTATAACGACAATACTGATTGTATATCTACGTCTTGGTCGTTTTCAAATACTTGTTTTCCGTTCAGTTCCTTTAACTCAATGATGAAGTTTACGTATACCTTTTTCGGATGGAGCTCTTTCACGAGGTTGCAGGCGGCTTTCATAGTACCTCCGGTAGCCAGTAAGTCATCGTGAGACAATACGACGTCGTTCTCGTTCAATGCATCTTTGTGTATCTGCACGGTGTCTTT        GCF_000025985.1_0

caclulating embeddings -> generates three files:
(base) andyfely@Andys:/mnt/c/Users/andyf/SWE/data/short_reads/embeddings$ ls
mock_test_emb.npy  mock_test_emb_stand.npy  mock_test_labels.txt

clsutering embeddings

visualize embendings

nextflow pipeline: from input.txt to visualization

new scripts

- pipeline.nf
- visualize_embeddings.py
- cluster_embeddings.py
- checkout_npy.py

modified scripts for pipeline:

- calculate_embedding_for_tsv.py
- utils_parallelized.py
- fastq_to_tsv_local.py
- download_genome.py
- process_genomes.py

### update: pull request changes von 31.01.2025 #1
modified scripts:
- pipeline.nf
- fastq_to_tsv.py

new files:
parallized_pipeline.nf
test_badread.sh

pipeline.nf:
included long read simulation process and adjusted processes if nessacry

fastq_to_tsv.py:
modified in order to correctly handle long read discriptors (diffrent formating that for short reads)

parallized_pipeline.nf
mostly the same as pipeline.nf but now parallized processes for short and long reads

test_badread.sh
had issues with badread tested diffrent parameters 

### update: pull request changes von 31.01.2025 #2
modified:
- parallized_pipeline.nf

new file:
- nextflow.config

-> created config file for netflow with the parameters and custome testing profile with parameter set for short runtime specifically

with testing parameters
usage: nextflow run parallized_pipeline.nf -profile test

usage: nextflow run parallized_pipeline.nf
