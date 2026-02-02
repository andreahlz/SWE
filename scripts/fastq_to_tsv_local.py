from Bio import SeqIO
from pathlib import Path
import os
import sys

def fastq_to_clustering_tsv(fastq_file):
    output_file = Path(fastq_file).stem.removesuffix("_genomic_R1")+'.tsv'
    with open(output_file,"w") as file:
        print("current tsv file in ",output_file)
        for record in SeqIO.parse(fastq_file,"fastq"):
            file.write(str(record.seq)+"\t"+record.id.split("_")[1]+"\n")


print("path to fastq file:",sys.argv[1])
fastq_to_clustering_tsv(sys.argv[1])
