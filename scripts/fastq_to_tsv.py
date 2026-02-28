"""
fastq_to_tsv.py
Converts fastq file into tsv with the colums: sequence, species. 
Output file name is given by command line argument.
Written for fastq files created by InSilicoSeq or Badread.

Usage
-----
    python fastq_to_tsv.py --input path/to/reads.fastq
    python fastq_to_tsv.py --input path/to/reads.fastq --output path/to/output.tsv

Arguments
---------
--input  : Path to input FASTQ file (required). Produced by ISS or Badread.
--output : Path for output TSV file (optional).
           Default: derived from input filename by stripping '_genomic_R1' suffix
           and replacing the extension with .tsv.


"""
from Bio import SeqIO
from pathlib import Path
import sys
import argparse
import re

# Regex to extract the accession ID from InSilicoSeq read IDs.
#iss hangs _number_number and /1 or /2 this pattern is not part of the accesion id and should not be included
iss_re = re.compile(r'^(.+?)_\d+_\d+(?:/\d+)?$')

def fastq_to_clustering_tsv(fastq_file,output_file=None):
    """
    Parses a FASTQ file and writes a two-column TSV (sequence, bin_id).
    Detects the read format automatically:
    - Badread (long reads): description contains whitespace.
    - InSilicoSeq (short reads): no whitespace, ID matched via iss_re.

    Paramaters:
    fastq_file : Path to the input FASTQ file.
    output_file: Path for the output
    """
    
    fastq_path = Path(fastq_file)
    if output_file is None:
        output_file = fastq_path.stem.removesuffix("_genomic_R1")+'.tsv'
    
    with open(output_file,"w") as file:
        #header of tsv
        file.write("sequence\tbin_id\n")
        #Note check how to handle the print statments
        print("current tsv file in ",output_file)

        for record in SeqIO.parse(str(fastq_path),"fastq"):
            # Extract label based on format
            # Short reads: @GCF_000025985.1_0_0/1
            # Long reads: @uuid GCF_000025985.1,+strand,start-end ..
            description = record.description #saves header
            if ' ' in description: #badread reads have whitespaces
                parts = description.split()
                if len(parts) >= 2:
                    genome_descriptor = parts[1]
                    label= genome_descriptor.split(',')[0]
                else:
                    raise ValueError(f"Unexpected Badread format: {description}")
            else: #in fastq produced by InSilicoSeq
                m = iss_re.match(record.id) #match pattern at the beginning of record.id
                if m:
                    label = m.group(1) #extract first captured group
                else: #no match was found
                     raise ValueError(f"Unexpected InSilicoSeq format: id={record.id} desc={description}")

            file.write(f"{record.seq}\t{label}\n")

    print(f"Written TSV: {output_file}", file=sys.stderr)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Convert FASTQ reads to TSV for embedding/clustering.")
    parser.add_argument("--input", required=True, help="Input FASTQ file")
    parser.add_argument("--output", default=None, help="Output TSV file (default: derived from input name)")
    args = parser.parse_args()

    fastq_to_clustering_tsv(args.input, args.output)
    