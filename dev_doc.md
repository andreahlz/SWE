# SWE Project Logbook

## NEXT TO DO:
- how to solve abbundance when multiple sequences
- covarage oder absolute zahl an reads oder was genau der unterschied
- long reads sim script

## NOTE neue eintrage über den alten (dann müssen wir nicht das ganze doc durchscrollen)

### update: pull request chnages von 13.12.2025

files changed:
- enviroment.yml
- simulate_short_reads.sh
- .gitignore

new files:
- config.env
- test_abbundance.txt 
- multifasta\test_mock_gut_genomes.fna

info on chnages:
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

