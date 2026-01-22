# SWE Project Logbook

## NEXT TO DO:
- how to solve abbundance when multiple sequences
- covarage oder absolute zahl an reads oder was genau der unterschied
- long reads sim script

## NOTE neue eintrage über den alten (dann müssen wir nicht das ganze doc durchscrollen)

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


### update: 21.1 and 22.1.
files changed: 
- config.env

new scripts:
- download_genomes.py
- process_genomes.py

deleted scripts:
- genomes_download.py

info on changes:
- config.env:
  added paragraphs: 
    - genome download parameters
    - genome processing parameters
    - adjusted genome_multifasta
    - adjusted abundance_file 
- download_genomes vs genomes_download
  new script 
    - doesnt have hardcoded paths -> config.env
    - keeps the output files separate -> doesnt combine them into multi fasta file
- process_genomes.py
  workflow:
    - Reads abundance values from your input.txt
    - For each .fna file in data/genomes/processed:
        - Extracts all headers
        - Combines multiple sequences into one long sequence
        - Creates new file with just the Accession ID as header
        - Saves to data/genomes/processed/newheader/
    - Creates a CSV file (header_info.csv) with:
        - Accession ID
        - Original headers (separated by semicolons if multiple)
        - Abundance value
    - Combines all processed files into one multi-fasta file
  output:
    - new folder: data/genomes/processed/newheader/
    - separate .fna files with new header (data/genomes/processed/newheader)
    - header_info_csv file (data/genomes)
    - combined_genomes.fna (data/genomes)



