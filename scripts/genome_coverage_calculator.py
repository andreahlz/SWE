#!/usr/bin/env python3
"""
Calculates the number of short reads required per genome to reach
params.target_coverage given the genome size and params.short_read_length.
Writes the total read count to coverage_short.txt.

coverage = (target_coverage * effektive_size)/short_read_length
effektive_size = len(sequence) * abundance

Usage: 
    python genome_coverage_calculator.py \
        --input_file            data/input/input.tx \
        --newheader_fastas_dir  processed \
        --target_coverage       10 \
        --short_read_length     300 \
        --txt_short_output      coverage_short.txt 

Arguments:
    --input_file            Tab-separated file with accession IDs and abundances (see Section 5)
    --newheader_fastas_dir  Directory containing the re-headered .fna files
    --target_coverage       How often each base is sequences on average
    --short_read_length     Length of short reads in bp 
    --txt_short_output      coverage_short.txt 
"""

#from dotenv import dotenv_values
from pathlib import Path
import argparse

def read_abundances(input_file):
    """Read abundance values from input.txt"""
    abundances = {}
    with open(input_file) as f:
        next(f)  
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                accession = parts[0]
                abundance = float(parts[1])
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
                
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))
            
                current_header = line[1:]  
                current_seq = []
            else:
                current_seq.append(line)
        
        
        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))
    
    return sequences


def get_genome_size(fna_file):
    """Get genome size by parsing FASTA file"""
    sequences = parse_fasta(fna_file)
    if not sequences:
        return 0
    # Return length of combined sequence 
    return len(sequences[0][1])


def main(args):
    
    """
    config = dotenv_values("../config.env")
    
    # Get parameters
    input_file = Path("..") / config["input_file"]
    newheader_dir = Path("..") / config["newheader_dir"]
    target_coverage = int(config["target_coverage"])
    short_read_length = int(config["short_read_length"])
    long_read_length = int(config["mean_length"]) 
    txt_short_output = Path("..") / config["coverage_short_txt"] 
    txt_long_output = Path("..") / config["coverage_long_txt"] 
    """
    input_file = Path(args.input_file)
    newheader_dir = Path(args.newheader_fastas_dir)
    target_coverage = int(args.target_coverage)
    short_read_length = int(args.short_read_length)
    txt_short_output = Path(args.txt_short_output) 
    abundances = read_abundances(input_file)
    
    genomes = {}
    for fna_file in sorted(newheader_dir.glob("*.fna")):
        accession = fna_file.stem
        size = get_genome_size(fna_file)
        abundance = abundances.get(accession, 0.0)
        genomes[accession] = {"size": size, "abundance": abundance}
    
    print("=" * 60)
    print(f"METAGENOMIC READ CALCULATION FOR {target_coverage}X COVERAGE")
    print("=" * 60)
    
    f_short = open(txt_short_output, 'w', newline='')
    #f_long = open(txt_long_output, 'w', newline='')

    print(f"\n### SHORT READS ({short_read_length} bp) ###")
    total_short_reads = 0
    for accession, data in genomes.items():
        effective_size = data["size"] * data["abundance"]
        reads_needed = (target_coverage * effective_size) / short_read_length
        total_short_reads += reads_needed

        print(f"{accession:40} {int(reads_needed):>10,} reads")
    f_short.write(f"{total_short_reads:.0f} \n")

    print(f"\n{'TOTAL SHORT READS':40} {int(total_short_reads):>10,} reads")

    """
     print(f"\n### LONG READS ({long_read_length/1000:.0f} kb) ###")
    total_long_reads = 0
    for accession, data in genomes.items():
        effective_size = data["size"] * data["abundance"]
        reads_needed = (target_coverage * effective_size) / long_read_length
        total_long_reads += reads_needed

        print(f"{accession:40} {int(reads_needed):>10} reads")
        f_long.write(f"{accession:40} {int(reads_needed):>10} \n")

    
    print(f"\n{'TOTAL LONG READS':40} {int(total_long_reads):>10} reads")
    """
   

  
    print(f"Saved coverage info for simulating short reads to {txt_short_output}")
    #print(f"Saved coverage info for simulating long reads to {txt_long_output}")

    f_short.close()
    #f_long.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate coverage based read numbers')
    parser.add_argument('--input_file', type=str, required=True,
                        help="Input file with accesions and abbundance")
    parser.add_argument('--newheader_fastas_dir',type=str, required=True,
                        help="directory with processesed fasta files (with adjusted header and combined sequences)")
    parser.add_argument('--target_coverage',type=int, required=True,
                        help="Desired sequencing coverage")
    parser.add_argument('--short_read_length',type=int, required=True,
                        help="length of short reads (MiSeq = 300bp, HiSeq = 125bp, NextSeq = 300bp, Novaseq = 150bp)")
    parser.add_argument('--txt_short_output',type=str,required=True,help='number of reads for correct coverage')
    args = parser.parse_args()
    main(args)