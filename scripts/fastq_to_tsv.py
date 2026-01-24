from Bio import SeqIO
from pathlib import Path
import sys
import argparse


def fastq_to_clustering_tsv(fastq_file,output_file=None):
    fastq_path = Path(fastq_file)
    if output_file is None:
        output_file = fastq_path.stem.removesuffix("_genomic_R1")+'.tsv'
    
    with open(output_file,"w") as file:
        #header of tsv
        file.write("sequence\tbin_id\n")
        #Note check how to handle the print statments
        print("current tsv file in ",output_file)

        for record in SeqIO.parse(str(fastq_path),"fastq"):
            parts = record.id.split("_")
            if len(parts) < 2:
                raise ValueError(
                    f"Read ID '{record.id}' does not contain '_' -> cannot extract label with split('_')[1]"
                )
            label = parts[1]
            file.write(f"{record.seq}\t{label}\n")

    print(f"Written TSV: {output_file}", file=sys.stderr)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Convert FASTQ reads to TSV for embedding/clustering.")
    parser.add_argument("--input", required=True, help="Input FASTQ file")
    parser.add_argument("--output", default=None, help="Output TSV file (default: derived from input name)")
    args = parser.parse_args()

    fastq_to_clustering_tsv(args.input, args.output)
    