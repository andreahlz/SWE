#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="${script_dir%/scripts}"

# Relative paths to project-root
genome_input_file="${project_root}/data/genomes/mock_gut.fna"
abundance_input_file="${project_root}/data/abundances/abundance.txt"
output_dir="${project_root}/data/simulated_reads/long_reads"
file_name="test_01"

# Badread parameters
mean_length=10000        # Mean read length (ONT typical)
length_stdev=5000        # Standard deviation
mean_identity=95         # Mean identity (95% = 5% error, ONT typical)
identity_stdev=5         # Standard deviation for identity
seed=42
cpus=8

mkdir -p "${output_dir}"
temp_dir="${output_dir}/temp"
mkdir -p "${temp_dir}"

# Read abundances
declare -A abundances
while IFS=$'\t' read -r species abundance; do
    [[ "$species" == "species" ]] && continue  # Skip header
    abundances["$species"]=$abundance
done < "${abundance_input_file}"

# Split multifasta into individual genomes
echo "Splitting multifasta into individual genomes..."
awk '/^>/ {if (seq) print seq; seq=""; split($0,a," "); gsub(/^>/,"",a[1]); filename=a[1]".fna"; print $0 > "'${temp_dir}'/"filename; next} {seq=seq$0} END {if (seq) print seq}' "${genome_input_file}"

# Simulate reads for each genome based on abundance
echo "Simulating long reads..."
for genome_file in "${temp_dir}"/*.fna; do
    [[ ! -f "$genome_file" ]] && continue
    
    genome_name=$(basename "$genome_file" .fna)
    abundance=${abundances[$genome_name]:-0.01}  # Default if not found
    
    # Calculate coverage based on abundance (scale to get reasonable read counts)
    coverage=$(echo "$abundance * 10" | bc)  # 10x base coverage
    
    echo "  Processing: $genome_name (abundance: $abundance, coverage: ${coverage}x)"
    
    badread simulate \
        --reference "${genome_file}" \
        --quantity "${coverage}x" \
        --length ${mean_length},${length_stdev} \
        --identity ${mean_identity},${identity_stdev},${identity_stdev} \
        --seed ${seed} \
        > "${temp_dir}/${genome_name}.fastq" 2>/dev/null
done

# Combine all reads into single output
echo "Combining reads..."
cat "${temp_dir}"/*.fastq > "${output_dir}/${file_name}.fastq"

# Compress
echo "Compressing..."
gzip -f "${output_dir}/${file_name}.fastq"

# Cleanup
rm -rf "${temp_dir}"

echo "✓ Long read simulation complete: ${output_dir}/${file_name}.fastq.gz"