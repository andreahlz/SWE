"""
cluster_embeddings.py
K-means clustering for combined embeddings.

Based on eval_clustering_classification_changed.py,
modified for the Nextflow pipeline.

Loads pre-computed DNABERT-S embeddings (.npy) and their true genome labels (.txt),
runs K-means clustering across 5 random seeds, and saves the best result plus
mean metrics across all seeds.

Usage
-----
    python cluster_embeddings.py \\
        --embedding_file path/to/<base>_emb_stand.npy \\
        --label_file     path/to/<base>_labels.txt \\
        --output_file    path/to/<base>_kmeans_results.csv

Arguments
---------
--embedding_file : Path to standardised embedding matrix (.npy), shape (N_reads x embedding_dim).
                   Use the _emb_stand.npy output from calculate_embedding_for_tsv.py.
--label_file     : Path to labels file (.txt) produced by calculate_embedding_for_tsv.py.
                   First line must be the header "true_label".
--output_file    : Path for the output CSV. A companion metrics .txt file is written
                   to the same directory automatically (suffix _metrics.txt).

Notes:
- Embeddings are L2-normalised before clustering
- K-means is run with 5 seeds (0-4) and n_init=10. The seed with the highest ARI
  is selected as the best result for the CSV output.
- Labels are mapped to integer IDs via sorted(set(labels)) to ensure a consistent
  and reproducible mapping across runs.
- num_clusters is derived directly from the number of unique labels
"""

import argparse
import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import normalize
import sklearn.metrics
from sklearn.cluster import KMeans


def load_embeddings_from_npy(npy_file,label_file):
    """
    Load embedding matrix and true labels
    parameters:
        npy_file :Path to .npy file containing the embedding matrix (N_reads x embedding_dim).
        label_file : Path to labels .txt file. First line is the header and is skipped.
    returns:
    embeddings : np.ndarray, shape (N_reads, embedding_dim)
    labels : True genome labels, one per read.
    
    """
    embeddings = np.load(npy_file)
    with open(label_file, 'r') as f:

        lines = f.read().strip().split('\n')
        labels = lines[1:] #header skip
    return embeddings, labels


def prepare_labels(labels):
    """
    convert labels to numeric values (for sklearn compatibility)
    Parameters:
        labels: True genome labels (e.g. NCBI accession IDs).
    returns:
        labels_numeric : np.ndarray of int, Integer-encoded labels aligned with the input list.  
        num_clusters: Number of unique species / clusters
        label2id: Mapping from original label string to integer ID
    """  
    label2id = {l: i for i, l in enumerate(sorted(set(labels)))}
    labels_numeric = np.array([label2id[l] for l in labels])
    num_clusters = len(label2id)

    return labels_numeric, num_clusters, label2id

def perform_kmeans(embeddings, labels, num_clusters, random_seeds=[0, 1, 2, 3, 4]):
    """
    K-meand clustering across multiple random seeds
        
    Embeddings are L2-normalised before clustering to approximate cosine similarity.
    Each seed produces an independent clustering; all results are returned for
    downstream selection (best ARI) and averaging.
    
    Parameters:
        embeddings: np.ndarray, shape (N_reads, embedding_dim)
        labels: Integer-encoded true labels from prepare_labels().
        num_clusters: Number of clusters (equals number of unique species)
        random_seeds: Seeds used for KMeans random_state. Default: [0, 1, 2, 3, 4]
    
    returns:
        results : list of dict
        One dict per seed with keys: seed, predicted_labels, purity,
        completeness, ari, nmi.
    """
    embedding_norm = normalize(embeddings)
    results = []
    for random_seed in random_seeds:
        ### perform k-means clustering
        kmeans = KMeans(n_clusters=num_clusters,
                        random_state=random_seed,
                        max_iter=1000,
                        init="k-means++",
                        n_init=10)
        kmeans.fit(embedding_norm)
        preds_clustering = kmeans.labels_
        homogeneity = sklearn.metrics.homogeneity_score(labels, preds_clustering)
        completeness = sklearn.metrics.completeness_score(labels, preds_clustering)
        ari = sklearn.metrics.adjusted_rand_score(labels, preds_clustering)
        nmi = sklearn.metrics.normalized_mutual_info_score(labels, preds_clustering)

        print(f"Kmeans homogeneity: {homogeneity} completeness: {completeness} ari: {ari} nmi: {nmi}")

        #diff to old implementation: reuslts saved in dic instead numpy array
        results.append({
            'seed':random_seed,
            'predicted_labels': preds_clustering,
            'homogeneity': homogeneity,
            'completeness' : completeness,
            'ari':ari,
            'nmi':nmi})
    return results

def calculate_mean_results(results):
    """
    Compute mean clustering metrics across all seeds.
    Parameters:
        results : list of dict, Output of perform_kmeans(); each dict contains homogeneity, completeness, ari, nmi.
    Returns:
        mean_results :Mean value for each metric across all seeds
    """
    metrics = ['homogeneity', 'completeness','ari','nmi']
    mean_results = {
        metric: np.mean([r[metric] for r in results])
        for metric in metrics
    }
    return mean_results

def save_results(results, original_labels, output_file):
    """
    Clustering results saved to csv and metrics summary as txt file
    The best seed is selected by highest ARI. The CSV contains true and predicted
    labels for every read. The metrics file contains best-seed and mean metrics.
    Parameters:
        results :Output of perform_kmeans()
        original_labels: Original string genome labels (not integer-encoded)
        output_file: Path for the output.
    """
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
        f.write(f"Best Result (Seed {best_result['seed']}:)\n")
        f.write(f"ARI: {best_result['ari']:.4f}\n")
        f.write(f"homogeneity: {best_result['homogeneity']:.4f}\n")
        f.write(f"Completeness: {best_result['completeness']:.4f}\n")
        f.write("Mean Results (across all seeds):\n")
        for metric in ['homogeneity', 'completeness', 'ari', 'nmi']:
            f.write(f"  {metric.capitalize():13s}: {mean_results[metric]:.4f}\n")
        

def main(args):
    # Load data
    embeddings, original_labels = load_embeddings_from_npy(args.embedding_file, args.label_file)
    #check if embedding number matches labels number (check for missing data whhich would make results unreliable)
    if embeddings.shape[0] != len(original_labels):
        raise ValueError(f"Mismatch: embeddings = {embeddings.shape[0]} labels = {len(original_labels)}")
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
                        help='Path to embeddings.npy')
    parser.add_argument('--label_file', type=str, required=True,
                        help='Path to labels.txt')
    parser.add_argument('--output_file', type=str, required=True,
                        help='Path to save clustering results')
    
    args = parser.parse_args()
    main(args)

