#!/bin/bash
set -euo pipefail


script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="${script_dir%/scripts}"

mkdir -p "${project_root}/data_test/test_long"

fasta_file="${project_root}/data/short_reads/genomes/fasta_files/GCF_000025985.1.fna"

# Test OHNE --identity (nutzt Badread default)
badread simulate \
    --reference "${fasta_file}" \
    --quantity 0.5x \
    --length 5000,2000 \
    --seed 42 \
    > test.fastq

echo ""
echo "Generierte FASTQ-Datei: test.fastq"
echo "Anzahl Reads: $(grep -c "^@" test.fastq)"
echo ""
echo "Erste 8 Zeilen:"
head -n 8 test.fastq