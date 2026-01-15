#!/bin/bash
set -euo pipefail


script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="${script_dir%/scripts}"

source "${project_root}/config.env"


mkdir -p "${project_root}/${out_dir_shortreads}"

iss generate\
    --genomes "${project_root}/${genome_multifasta}" \
    --abundance_file "${project_root}/${abundance_file}" \
    --model "${error_model}" \
    --n_reads "${n_reads}" \
    --seed "${seed}" \
    --cpus "${cpus}" \
    --output "${project_root}/${out_dir_shortreads}/${short_reads_output}"