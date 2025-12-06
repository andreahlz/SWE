#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="${script_dir%/scripts}"

#relative paths to project-root
genome_input_file="${project_root}/data/genomes/mock_gut.fna" #Multifasta file dann
abundance_input_file="${project_root}/data/abundances/abundance.txt"
output_dir="${project_root}/data/simulated_reads/short_reads"
file_name="test_01"

#insilco parameters
n_reads=2400000
error_model="MiSeq"
seed=42
cpus=8

mkdir -p "${output_dir}"

iss generate\
    --genomes ${genome_input_file}\
    --abundance_file ${abundance_input_file}\
    --model ${error_model}\
    --n_reads ${n_reads}\
    --seed ${seed}\
    --cpus ${cpus}\
    --output ${output_dir}/${file_name}