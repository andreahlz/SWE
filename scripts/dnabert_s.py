from xml.parsers.expat import model
import torch
from transformers import AutoTokenizer, AutoModel
from pathlib import Path
import os
from Bio import SeqIO

###local test: file has another name
conf_filename = 'config.env'
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

def parse_reads_file(path_sequence_file):
    sequences = []
    for record in SeqIO.parse(path_sequence_file, "fastq"):    
        print(record)            
        sequences.append(str(record.seq))   
    print(f"DEBUG: Inside function, sequences has {len(sequences)} items")
    print(f"DEBUG: First sequence: {sequences[0] if sequences else 'EMPTY'}")
    return sequences

def calculate_embedding(sequences):
    tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-S", trust_remote_code=True)
    model = AutoModel.from_pretrained("zhihan1996/DNABERT-S", trust_remote_code=True)
    dna = "ACGTAGCATCGGATCTATCTATCGACACTTGGTTATCGATCTACGAGCATCTCGTTAGC"
    inputs = tokenizer(dna, return_tensors = 'pt')["input_ids"]
    hidden_states = model(inputs)[0] # [1, sequence_length, 768]
    # embedding with mean pooling
    embedding_mean = torch.mean(hidden_states[0], dim=0)
    print(embedding_mean.shape) # expect to be 768
    embeddings = []
    for seq in sequences:
        tok = tokenizer(seq, return_tensors="pt")["input_ids"]
        inputs = tokenizer(seq, return_tensors="pt")["input_ids"]
        hidden_states = model(inputs)[0]
        mod = model(inputs)[0]
        seq_embedding = torch.mean(hidden_states[0], dim=0)
        embeddings.append(seq_embedding)
    return embeddings
#######Muss ich noch ändern
proj_root = Path(__file__).parent.parent
conf_path = proj_root/conf_filename
config = load_config(conf_path)
#Short reads
short_read_file_path = proj_root/config['out_dir_shortreads']/config['short_reads_output']
short_read_sequences = parse_reads_file(short_read_file_path)
'''out_dir_shortreads="data/simulated_reads/short_reads"
out_dir_longreads="data/simulated_reads/long_reads"'''
####short_reads_output="short_reads"
##long_reads_output="long_reads"
#parse_reads_file()    
calculate_embedding(short_read_sequences)
print("test")  