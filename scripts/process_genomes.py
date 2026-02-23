#!/usr/bin/env python3
"""
This script processes genome files:
1. Extract headers and save to CSV
2. Combine sequences within each file
3. Add new headers with Accession ID
4. Create combined multi-fasta file

Usage: python process_genomes.py

Requirements: pip install python-dotenv
"""

import csv
import argparse
from pathlib import Path
#from dotenv import dotenv_values


def read_abundances(input_file):
    """Read abundance values from input.txt"""
    abundances = {}
    with open(input_file) as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                accession = parts[0]
                abundance = parts[1]
                abundances[accession] = abundance
    return abundances


def parse_fasta(fasta_file):
    """Parse a FASTA file and return list of (header, sequence) tuples"""
    sequences = []
    current_header = None
    current_seq = []
    
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))
                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))
    
    return sequences


def process_genome_file(fna_file, abundances, output_dir):
    """
    Process a single genome file:
    - Extract headers
    - Combine sequences
    - Save with new header
    """
    # Get accession ID from filename (e.g., GCF_000025985.1.fna -> GCF_000025985.1)
    accession = fna_file.stem
    
    # Parse the FASTA file
    sequences = parse_fasta(fna_file)
    
    if not sequences:
        print(f"  Warning: No sequences found in {fna_file.name}")
        return None
    
    # Extract all headers
    headers = [seq[0] for seq in sequences]
    
    # Combine all sequences
    combined_sequence = ''.join([seq[1] for seq in sequences])
    
    # Get abundance
    abundance = abundances.get(accession, "N/A")
    
    # Create new file with accession ID as header
    output_file = output_dir / f"{accession}.fna"
    with open(output_file, 'w') as f:
        f.write(f">{accession}\n")
        # Write sequence in lines of 80 characters (standard FASTA format)
        for i in range(0, len(combined_sequence), 80):
            f.write(combined_sequence[i:i+80] + '\n')
    
    print(f"  ✓ Processed {accession}: {len(sequences)} sequence(s) combined")
    
    return {
        'accession': accession,
        'headers': '; '.join(headers),  # Join multiple headers with semicolon
        'abundance': abundance
    }


def create_multifasta(input_dir, output_file):
    """Combine all .fna files into one multi-fasta file"""
    with open(output_file, 'w') as out:
        for fna_file in sorted(input_dir.glob("*.fna")):
            with open(fna_file) as f:
                out.write(f.read())
    print(f"✓ Created multi-fasta: {output_file}")


def main(args):
    '''
    # Load configuration
    config = dotenv_values("../config.env")
    
    # Get parameters
    input_file = Path("..") / config["input_file"]
    single_fastas_dir = Path("..") / config["processed_genomes_dir"]
    output_dir = Path("..") / config["newheader_dir"]
    csv_output = Path("..") / config["headers_csv"]
    multifasta_output = Path("..") / config["combined_multifasta"]
    '''
    input_file = Path(args.input_file)
    single_fastas_dir = Path(args.single_fastas_dir)
    output_dir = Path(args.output_dir)
    csv_output =Path(args.csv_output)
    multifasta_output = Path(args.multifasta_output)
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Step 1: Reading abundances...")
    abundances = read_abundances(input_file)
    print(f"  Found {len(abundances)} accessions with abundance values\n")
    
    print("Step 2: Processing genome files...")
    csv_data = []
    
    for fna_file in sorted(single_fastas_dir.glob("*.fna")):
        result = process_genome_file(fna_file, abundances, output_dir)
        if result:
            csv_data.append(result)
    
    print(f"\nStep 3: Creating CSV with header information...")
    csv_output.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['accession', 'headers', 'abundance'])
        writer.writeheader()
        writer.writerows(csv_data)
    print(f"  ✓ Saved header info to {csv_output}")
    
    print(f"\nStep 4: Creating combined multi-fasta file...")
    multifasta_output.parent.mkdir(parents=True, exist_ok=True)
    create_multifasta(output_dir, multifasta_output)
    
    print(f"\n✓ All done! Processed {len(csv_data)} genome files")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Genome processing: from fasta to multifasta')
    parser.add_argument('--input_file', type=str, required=True,
                        help="Input file with accesions and abbundance")
    parser.add_argument('--single_fastas_dir', type=str, required=True,
                        help='directory were the fasta files are')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='output dir for fasta with new header')
    parser.add_argument('--csv_output', type=str, required=True,
                        help='path for csv_output')
    parser.add_argument('--multifasta_output', type=str, required=True,
                        help='name/path for multifasta_output')
    
    args = parser.parse_args()
    main(args)
