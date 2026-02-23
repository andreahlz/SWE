#!/usr/bin/env python3
"""
This script downloada genomes from NCBI by accession number
and reads accession IDs from input.txt

Usage: python genome_download.py

Requirements: pip install requests python-dotenv
"""

import sys
import zipfile
import io
from pathlib import Path
import argparse

try:
    import requests
except ImportError:
    print("ERROR: requests library not found")
    print("Install with: pip install requests")
    sys.exit(1)

try:
    from dotenv import dotenv_values
    HAS_DOTENV = True
except ImportError:
    HAS_DOTENV = False
    print("ERROR: python-dotenv library not found")
    print("Install with: pip install python-dotenv")
    sys.exit(1)


def download_genome(accession, output_dir, timeout):
    """Download a single genome by accession using NCBI API"""
    print(f"Downloading {accession}...")
    
    # NCBI datasets API URL
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"
    params = {"include_annotation_type": "GENOME_FASTA"}
    
    try:
        # Download the zip file
        response = requests.get(url, params=params, stream=True, timeout=timeout)
        response.raise_for_status()
        
        # Read zip from memory
        zip_data = io.BytesIO(response.content)
        
        # Extract FASTA files
        with zipfile.ZipFile(zip_data) as zf:
            # Find .fna files
            fna_files = [f for f in zf.namelist() if f.endswith('.fna')]
            
            if not fna_files:
                print(f"No FASTA file found in download")
                return False
            
            # Extract first .fna file
            fna_content = zf.read(fna_files[0])
            
            # Save to output directory
            output_fasta = output_dir / f"{accession}.fna"
            output_fasta.write_bytes(fna_content)
            
            print(f"Saved to {output_fasta}")
            return True
        
    except requests.exceptions.RequestException as e:
        print(f"Download failed: {e}")
        return False
    except Exception as e:
        print(f"Failed: {e}")
        return False


def main(args):
    if args.input_file is None and HAS_DOTENV:
    # Get parameters from config (adjust paths to be relative to project root)
        config = dotenv_values("../config.env")
        input_file = Path("..") / config["input_file"]
        output_dir = Path("..") / config["genome_output_dir"]
        timeout = int(config["download_timeout"])
    else:
        input_file = Path(args.input_file)
        output_dir = Path(args.output_dir)
        timeout = int(args.timeout)
 
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read accession IDs from input.txt
    accessions = []
    with open(input_file) as f:
        # Skip header line
        next(f)
        for line in f:
            line = line.strip()
            if line:
                # Extract first column (Accession_ID)
                accession = line.split()[0]
                accessions.append(accession)
    
    print(f"Found {len(accessions)} accessions in {input_file}")
    print(f"Output directory: {output_dir}")
    
    # Download each genome
    success = 0
    for acc in accessions:
        if download_genome(acc, output_dir, timeout):
            success += 1
    
    print(f"\nDownloaded {success}/{len(accessions)} genomes")
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download gemones from NCBI')
    parser.add_argument('--input_file', required=True, help="Input file with accesions and abbundance")
    parser.add_argument('--output_dir', required=True, help="output directory for single fasta files")
    parser.add_argument('--timeout', type=int, default=300, help='Download timeout in seconds')
    args = parser.parse_args()
    main(args)