def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

import argparse
import os
from sklearn.preprocessing import normalize
import csv
import sys
import numpy as np
import sklearn.metrics
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from pathlib import Path
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
import pandas as pd

from utils_parallelized import get_embedding

DEBUG_MODE = False

def main(args):
    model_list = args.model_list.split(",")
    for model in model_list:
        #output csv should include: embedding standard,embedding,labels
        # -> have to make 2 files
        embeddings_all = {}
        labels_data = {}
        colors_datapoints = {}
        tsv_file_path = args.tsv_file_path
        tsv_file=os.path.basename(tsv_file_path)
        print(f"Start {model}  clustering")
        data_file = os.path.join(args.data_dir, f"tsvs/{tsv_file}")
                        
        with open(tsv_file_path, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            data = list(reader)[1:]

        dna_sequences = [d[0] for d in data]
        labels = [d[1] for d in data]
        # generate embedding
        embedding = get_embedding(dna_sequences, model,tsv_file, test_model_dir=args.test_model_dir)
        embedding_norm = normalize(embedding)
        #use embedding_standard for t-SNE
        embedding_standard = StandardScaler().fit_transform(embedding)
        #converts the values of the embedding into differences to the mean in standard deviations
        #reduce dimensionality of current embedding, save low dimensional embedding in dictionary
        '''df_embedding = pd.DataFrame({
                'embedding':list(embedding),
                'labels':labels
            })
        df_emb_stand=pd.DataFrame({
            'standard':list(embedding_standard),
            'labels':labels
        })
        df_embedding.to_csv(f"data/csv_embeddings/{Path(tsv_file).stem}.csv",header=False,index=False)
        df_emb_stand.to_csv(f"data/csv_embeddings/{Path(tsv_file).stem}_stand.csv",header=False,index=False)
        print("test")'''
        #instead of csv using npy

        output_dir = "data/species_embeddings"
        os.makedirs(output_dir, exist_ok=True)
        base_name = Path(tsv_file).stem
        np.save(f"{output_dir}/{base_name}_emb.npy",  embedding)
        np.save(f"{output_dir}/{base_name}_emb_stand.npy", embedding_standard)
        with open(f"{output_dir}/{base_name}_labels.txt", 'w') as f:
            f.write('\n'.join(labels))

if __name__ == "__main__":
    os.environ["TOKENIZERS_PARALLELISM"] = "false"
    if DEBUG_MODE:
        #muss ich noch ausbauen
        proj_root = Path(__file__).parent.parent
        args = argparse.Namespace(
            #todo: in bash skript download von DNABERTS starten, Pfad bei test_model_dir übergeben
            test_model_dir = "/home/barbara/evaluate/DNABERT-S",
            tsv_file = "clustering_short.tsv",
            model_list="test",
            data_dir="data"       
        )
    else:
        parser = argparse.ArgumentParser(description='Evaluate clustering')
        parser.add_argument('--model_list', type=str, default="tnf, test", help='List of models to evaluate, separated by comma. Currently support [tnf, tnf-k, dnabert2, hyenadna, nt, test]')
        parser.add_argument('--data_dir', type=str, default="/root/data", help='Data directory')
        # added sample size
        parser.add_argument('--sample_number', type=int, default=1, help='Number of samples to evaluate')
        parser.add_argument('--tsv_file_path',type = str, default="test.tsv",help='Path to tsv file')
        parser.add_argument('--test_model_dir',type=str,default='/root/evaluate/DNABERT-S',help='Directory to save trained models to test')
        args = parser.parse_args()
    print(args)
    print(type(args))
    main(args) 
