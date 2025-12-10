#!/usr/bin/env python3
"""
Download genomes from NCBI by accession number
Usage: python genome_download.py accessions.txt

Requirements: pip install requests
"""

import sys
import zipfile
import io
from pathlib import Path

try:
    import requests
except ImportError:
    print("ERROR: requests library not found")
    print("Install with: pip install requests")
    sys.exit(1)

def download_genome(accession, output_dir):
    """Download a single genome by accession using NCBI API"""
    print(f"Downloading {accession}...")
    
    # NCBI datasets API URL
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"
    params = {"include_annotation_type": "GENOME_FASTA"}
    
    try:
        # Download the zip file
        response = requests.get(url, params=params, stream=True, timeout=300)
        response.raise_for_status()
        
        # Read zip from memory
        zip_data = io.BytesIO(response.content)
        
        # Extract FASTA files
        with zipfile.ZipFile(zip_data) as zf:
            # Find .fna files
            fna_files = [f for f in zf.namelist() if f.endswith('.fna')]
            
            if not fna_files:
                print(f"  ✗ No FASTA file found in download")
                return False
            
            # Extract first .fna file
            fna_content = zf.read(fna_files[0])
            
            # Save to output directory
            output_fasta = output_dir / f"{accession}.fna"
            output_fasta.write_bytes(fna_content)
            
            print(f"  ✓ Saved to {output_fasta}")
            return True
        
    except requests.exceptions.RequestException as e:
        print(f"  ✗ Download failed: {e}")
        return False
    except Exception as e:
        print(f"  ✗ Failed: {e}")
        return False

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 download_genomes.py accessions.txt")
        print("\nFormat of accessions.txt:")
        print("GCF_000005845.2  # E. coli")
        print("GCF_000009045.1  # B. subtilis")
        sys.exit(1)
    
    accession_file = Path(sys.argv[1])
    output_dir = Path("data/genomes/processed")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read accessions
    accessions = []
    with open(accession_file) as f:
        for line in f:
            line = line.split('#')[0].strip()  # Remove comments
            if line:
                accessions.append(line)
    
    print(f"Found {len(accessions)} accessions")
    print(f"Output: {output_dir}\n")
    
    # Download each genome
    success = 0
    for acc in accessions:
        if download_genome(acc, output_dir):
            success += 1
    
    print(f"\n✓ Downloaded {success}/{len(accessions)} genomes")
    
    # Create combined FASTA
    combined = output_dir.parent / "combined_genomes.fna"
    with open(combined, 'w') as out:
        for fasta in output_dir.glob("*.fna"):
            out.write(open(fasta).read())
    
    print(f"✓ Combined FASTA: {combined}")

if __name__ == "__main__":
    main()