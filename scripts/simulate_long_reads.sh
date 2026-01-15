#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="${script_dir%/scripts}"

# Source config file
source "${project_root}/config.env"

# Set paths using config variables
genome_input_file="${project_root}/${genome_multifasta}"
abundance_input_file="${project_root}/${abundance_file}"
output_dir="${project_root}/${out_dir_longreads}"
file_name="${long_reads_output}"

# Create output directory
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
awk '/^>/ {
    if (seq) print seq; 
    seq=""; 
    split($0,a," "); 
    gsub(/^>/,"",a[1]); 
    filename=a[1]".fna"; 
    print $0 > "'${temp_dir}'/"filename; 
    next
} 
{seq=seq$0} 
END {if (seq) print seq}' "${genome_input_file}"

# Simulate reads for each genome based on abundance
echo "Simulating long reads with Badread..."
for genome_file in "${temp_dir}"/*.fna; do
    [[ ! -f "$genome_file" ]] && continue
    
    genome_name=$(basename "$genome_file" .fna)
    abundance=${abundances[$genome_name]:-0.01}
    
    # Calculate coverage based on abundance
    coverage=$(echo "$abundance * 10" | bc)
    
    echo "  ${genome_name}: abundance=${abundance}, coverage=${coverage}x"
    
    badread simulate \
        --reference "${genome_file}" \
        --quantity "${coverage}x" \
        --length ${mean_length},${length_stdev} \
        --identity ${mean_identity},${identity_stdev},${identity_stdev} \
        --seed ${seed} \
        > "${temp_dir}/${genome_name}.fastq" 2>/dev/null
done

# Combine all reads
echo "Combining reads..."
cat "${temp_dir}"/*.fastq > "${output_dir}/${file_name}.fastq"

# Compress
echo "Compressing..."
gzip -f "${output_dir}/${file_name}.fastq"

# Cleanup
rm -rf "${temp_dir}"

echo "✓ Long read simulation complete: ${output_dir}/${file_name}.fastq.gz"
