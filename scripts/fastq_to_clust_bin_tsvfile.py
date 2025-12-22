from Bio import SeqIO
from pathlib import Path

conf_filename = "config_dnaberts_test.env"

def fastq_to_clustering_tsv(fastq_file, output_tsv):
    path = Path(proj_root/f"data/DNABERTS_input/clustering_{output_tsv}.tsv")
    with path.open("w") as file:
        file.write("sequence\tbin_id\n")
        for record in SeqIO.parse(fastq_file,"fastq"):
            file.write(str(record.seq)+"\t"+record.id.split("_")[0]+"\n")

def load_config(config_path):
    """Lade Config-File und gib Dictionary zurück"""
    config = {}
    with open(config_path, 'r') as f:
        for line in f:
            line = line.strip()
            # Überspringe Kommentare und leere Zeilen
            if line.startswith('#') or not line or '=' not in line:
                continue 
            # Parse key=value
            key, value = line.split('=', 1)
            # Entferne Anführungszeichen
            value = value.strip('"').strip("'")
            config[key] = value    
    return config

proj_root = Path(__file__).parent.parent
conf_path = proj_root/conf_filename
configs = load_config(conf_path)
output_fastqs = [[configs["short_reads_output"],configs["out_dir_shortreads"],"short"],[configs["long_reads_output"],configs["out_dir_longreads"],"long"]]
path_output_dir = Path(proj_root/"data/DNABERTS_input")
path_output_dir.mkdir(exist_ok=True,parents=True)
for output_fastq_file in output_fastqs:
    path_output_fastq = proj_root/output_fastq_file[1]/output_fastq_file[0]
    fastq_to_clustering_tsv(path_output_fastq,output_fastq_file[2])
print("test")
        



    




















































"""
    Convert a FASTQ file to a clustering TSV file format.

    Parameters:
    fastq_file (str): Path to the input FASTQ file.
    output_tsv (str): Path to the output TSV file.
    
    with open(output_tsv, 'w') as tsv_out:
        # Write header
        tsv_out.write("Sequence_ID\tSequence\tQuality_Score\n")
        
        for record in SeqIO.parse(fastq_file, "fastq"):
            seq_id = record.id
            sequence = str(record.seq)
            quality_scores = ' '.join(map(str, record.letter_annotations["phred_quality"]))
            
            # Write to TSV
            tsv_out.write(f"{seq_id}\t{sequence}\t{quality_scores}\n")"""