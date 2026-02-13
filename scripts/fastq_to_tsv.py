from Bio import SeqIO
from pathlib import Path
import sys
import argparse
import re

#iss hangs _number_number and /1 or /2 this pattern is not part of the accesion id and should not be included
iss_re = re.compile(r'^(.+?)_\d+_\d+(?:/\d+)?$')

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
            # Extract label based on format
            # Short reads: @GCF_000025985.1_0_0/1
            # Long reads: @uuid GCF_000025985.1,+strand,start-end ..
            """"
            parts = record.id.split("_")
            if len(parts) < 2:
                raise ValueError(
                    f"Read ID '{record.id}' does not contain '_' -> cannot extract label with split('_')[1]"
                )
            label = parts[1]
            """
            description = record.description #saves header
            if ' ' in description: #badread reads have whitespaces
                parts = description.split()
                if len(parts) >= 2:
                    genome_descriptor = parts[1]
                    label= genome_descriptor.split(',')[0]
                else:
                    raise ValueError(f"Unexpected Badread format: {description}")
            else:
                m = iss_re.match(record.id)
                if m:
                    label = m.group(1)
                else:
                     raise ValueError(f"Unexpected InSilicoSeq format: id={record.id} desc={description}")

            file.write(f"{record.seq}\t{label}\n")

    print(f"Written TSV: {output_file}", file=sys.stderr)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Convert FASTQ reads to TSV for embedding/clustering.")
    parser.add_argument("--input", required=True, help="Input FASTQ file")
    parser.add_argument("--output", default=None, help="Output TSV file (default: derived from input name)")
    args = parser.parse_args()

    fastq_to_clustering_tsv(args.input, args.output)
    