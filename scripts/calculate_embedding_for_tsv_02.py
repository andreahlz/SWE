"""
Calculate embeddings for a TSV (sequence + label) and save:
- <base>_emb.npy          (raw embeddings)
- <base>_emb_stand.npy    (standardized embeddings)
- <base>_labels.txt       (labels, one per line)

Fixed for Nextflow:
- writes into a user-provided --out_dir (default: current directory)
- robust model_list parsing (trims whitespace)
- safer TSV parsing (skips header, ignores empty/broken rows)
"""

import argparse
import os
import csv
import numpy as np
from sklearn.preprocessing import StandardScaler

from utils_parallelized import get_embedding


def main():
    
    parser = argparse.ArgumentParser(description='Evaluate clustering')
    parser.add_argument('--model_list', type=str, default="tnf, test", help='List of models to evaluate, separated by comma. Currently support [tnf, tnf-k, dnabert2, hyenadna, nt, test]')
    parser.add_argument('--data_dir', type=str, default="/root/data", help='Data directory')   
    parser.add_argument('--tsv_file_path',type = str, default="test.tsv",help='Path to tsv file')
    parser.add_argument('--test_model_dir',type=str,default='/root/evaluate/DNABERT-S',help='Directory to save trained models to test')
    parser.add_argument('--out_dir', default='.', help='Output directory (default: current directory)')
    args = parser.parse_args()

    model_list = [m.strip() for m in args.model_list.split(",") if m.strip()]
    
    tsv_file_path = args.tsv_file_path
    tsv_file = os.path.basename(tsv_file_path)
    base_name = os.path.splitext(tsv_file)[0]

    
    output_dir = args.out_dir
    os.makedirs(output_dir, exist_ok=True)

    dna_sequences = []
    labels = []
    
    with open(tsv_file_path, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise ValueError(f"TSV is empty:{tsv_file_path}")
        for row in reader:
            if not row or len(row) <2:#scipng rows if invalid
                continue
            dna_sequences.append(row[0].strip())#saves sequence
            labels.append(row[1].strip())#saves labels
    model = model_list[0]
    # generate embedding
    embedding = get_embedding(dna_sequences, model,tsv_file, test_model_dir=args.test_model_dir, path_data_dir=args.data_dir)
    #use embedding_standard for t-SNE
    embedding_standard = StandardScaler().fit_transform(embedding)
    #converts the values of the embedding into differences to the mean in standard deviations
    #reduce dimensionality of current embedding, save low dimensional embedding in dictionary

    #instead of csv using npy
    np.save(os.path.join(output_dir, f"{base_name}_emb.npy"),  embedding)
    np.save(os.path.join(output_dir, f"{base_name}_emb_stand.npy"), embedding_standard)
    labels_path = os.path.join(output_dir, f"{base_name}_labels.txt")
    with open(labels_path, 'w') as f:
        for label in labels:
            f.write(f"{label}\n")

if __name__ == "__main__":
    main()
