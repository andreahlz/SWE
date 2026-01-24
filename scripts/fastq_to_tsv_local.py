from Bio import SeqIO
from pathlib import Path
import os
import sys

conf_filename = "config_dnaberts_test.env"

# def fastq_to_clustering_tsv(fastq_file):
#     output_file = Path(fastq_file).stem.removesuffix("_genomic_R1")+'.tsv'
#     tsv_dir = Path(Path(__file__).parent.parent/"data/tsvs")
#     tsv_dir.mkdir(parents=True,exist_ok=True)
#     with open(Path(tsv_dir/output_file),"w") as file:
#         file.write("sequence\tbin_id\n")
#         for record in SeqIO.parse(fastq_file,"fastq"):
#             file.write(str(record.seq)+"\t"+record.id.split("_")[1]+"\n")

# print("path to fastq file:",Path(__file__).parent.parent/"data/fastq_files"/sys.argv[1])
# fastq_to_clustering_tsv(Path(__file__).parent.parent/"data/fastq_files"/sys.argv[1])

# def fastq_to_clustering_tsv(fastq_file):
#     # output_file = Path(fastq_file).stem.removesuffix("_genomic_R1")+'.tsv'
#     output_file = Path(fastq_file).stem+'.tsv'
#     tsv_dir = Path(fastq_file).parent.parent/"tsvs"
#     print(tsv_dir)
#     tsv_dir.mkdir(parents=True,exist_ok=True)
#     with open(Path(tsv_dir/output_file),"w") as file:
#         print("current tsv file in ",Path(tsv_dir/output_file))
#         file.write("sequence\tbin_id\n")
#         for record in SeqIO.parse(fastq_file,"fastq"):
#             file.write(str(record.seq)+"\t"+record.id.split("_")[1]+"\n")

def fastq_to_clustering_tsv(fastq_file):
    output_file = Path(fastq_file).stem.removesuffix("_genomic_R1")+'.tsv'
    #output_file = Path(fastq_file).stem+'.tsv'
    with open(output_file,"w") as file:
        print("current tsv file in ",output_file)
        for record in SeqIO.parse(fastq_file,"fastq"):
            file.write(str(record.seq)+"\t"+record.id.split("_")[1]+"\n")


print("path to fastq file:",sys.argv[1])
fastq_to_clustering_tsv(sys.argv[1])