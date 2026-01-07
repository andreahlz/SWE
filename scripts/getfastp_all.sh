#!/bin/bash

fasta_dir="/home/barbara/SWE/data/genomes/processed"
fastp_dir="/home/barbara/SWE_trash/fastp_files"
for file in "$fasta_dir"/*; do
    echo "$file"
    filename="${file##*/}"
    test_output="/home/barbara/SWE_trash/fastp_files/$filename"
    iss generate --genomes "$file" --model MiSeq --n_reads 200000 --output "/home/barbara/SWE_trash/fastp_files/$filename"
done

