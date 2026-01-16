# SWE Project Logbook

## NEXT TO DO:
- testing pipline
- extend pipline from fasta to results 
- handling of embeddings currently emmebding ["string"],"lable" -> current solutuion in combined_embeddings.py and vizualisation.py -> look into better solution 

## NOTE neue eintrage über den alten (dann müssen wir nicht das ganze doc durchscrollen)

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
  ↓ (parallel!)
Calculate embeddings (parallel per species)
  ↓
Combine embeddings 
  ↓
K-means clustering 
  ↓
t-SNE visualization 
  ↓
PNG output + CSV results + metrics file

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




