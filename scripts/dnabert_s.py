import torch
from transformers import AutoTokenizer, AutoModel
from pathlib import Path
import os

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
    with open(path_sequence_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith('>'):
                sequences.append(line)
    return sequences
def calculate_embedding(path_sequence_file):
    print("test")  

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
print("test")  