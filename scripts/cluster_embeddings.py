"""
K-means clustering for combined embeddings

based on eval_clustering_classification_changed.py
-> modified for nextflow pipeline
"""

import argparse
import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import normalize
import sklearn.metrics
from sklearn.cluster import KMeans


def load_embeddings_from_csv(csv_file):

    '''
    load combined embeddings (CSV)
    Note: how to best handle embeddings?

    INPUT FORMAT NOW:
        Column 0: embedding as string "[0.1, 0.2, ...]"
        Column 1: label (species name)
    '''
    df = pd.read_csv(csv_file, header = None)
    embeddings_str = df[0].values 
    labels = df[1].values

    embeddings = []
    for emb_str in embeddings_str:
        #remove brackets nad parse to numpy array
        emb_str = emb_str.strip('[]')
        emb = np.fromstring(emb_str, sep= ' ')
        embeddings.append(emb)

    embeddings = np.array(embeddings)

    return embeddings,labels

def prepare_labels(labels):
    # convert labels to numeric values  
    label2id = {l: i for i, l in enumerate(set(labels))}
    labels_numeric = np.array([label2id[l] for l in labels])
    num_clusters = len(label2id)

    return labels_numeric, num_clusters, label2id

def perform_kmeans(embeddings, labels, num_clusters, random_seeds=[0, 1, 2, 3, 4]):
    #random_seeds = [0,1] #for quick testing, delete later
    embedding_norm = normalize(embeddings)
    results = []
    for random_seed in random_seeds:
        ### perform k-means clustering
        kmeans = KMeans(n_clusters=num_clusters,
                        random_state=random_seed,
                        max_iter=1000,
                        init="random",
                        n_init=3)
        kmeans.fit(embedding_norm)
        preds_clustering = kmeans.labels_
        purity = sklearn.metrics.homogeneity_score(labels, preds_clustering)
        completeness = sklearn.metrics.completeness_score(labels, preds_clustering)
        ari = sklearn.metrics.adjusted_rand_score(labels, preds_clustering)
        nmi = sklearn.metrics.normalized_mutual_info_score(labels, preds_clustering)

        print(f"Kmeans purity: {purity} completeness: {completeness} ari: {ari} nmi: {nmi}")

        #diff to old implementation: reuslts saved in dic instead numpy array
        results.append({
            'seed':random_seed,
            'predicted_labels': preds_clustering,
            'purity': purity,
            'completeness' : completeness,
            'ari':ari,
            'nmi':nmi})
    return results

def calculate_mean_results(results):
    metrics = ['purity', 'completeness','ari','nmi']
    mean_results = {
        metric: np.mean(r[metric] for r in results)
        for metric in metrics
    }
    return mean_results

def save_results(results, original_labels, output_file):
    best_result = max(results, key=lambda x: x['ari'])

    #CSV with results
    output_df = pd.DataFrame({
        'true_label':original_labels,
        'predicted_labels': best_result['predicted_labels']
    })

    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    output_df.to_csv(output_file, index=False)

    #report in txt
    metrics_file = output_file.replace('.csv', '_metrics.txt')
    mean_results = calculate_mean_results(results)
    with open(metrics_file,'w') as f:
        f.write("K-means Clustering Results\n")
        f.write(f"Best Result (Seed {best_result['seed']}:\n)")
        f.write(f"ARI: {best_result['ari']:.4f}\n")
        f.write(f"Purity: {best_result['purity']:.4f}\n")
        f.write(f"Completeness: {best_result['completeness']:.4f}\n")
        f.write(f"purity: {best_result['purity']:.4f}\n")
        f.write("Mean Results (across all seeds):\n")
        for metric in ['purity', 'completeness', 'ari', 'nmi']:
            f.write(f"  {metric.capitalize():13s}: {mean_results[metric]:.4f}\n")
        

def main(args):
    # Load data
    embeddings, original_labels = load_embeddings_from_csv(args.embedding_file)
    # Prepare labels
    labels_numeric, num_clusters,_ = prepare_labels(original_labels)
    # Clustering
    results = perform_kmeans(embeddings, labels_numeric, num_clusters)
    # Show mean
    mean_results = calculate_mean_results(results)
    # Save
    save_results(results, original_labels, args.output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='K-means clustering on combined embeddings')
    parser.add_argument('--embedding_file', type=str, required=True,
                        help='Path to combined_embeddings.csv')
    parser.add_argument('--output_file', type=str, required=True,
                        help='Path to save clustering results')
    
    args = parser.parse_args()
    main(args)

