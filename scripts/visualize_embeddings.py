"""
visualize_embeddings.py
t-SNE visualization of embeddings with clustering results

Based on: eval_clustering_classification_changed.py 
Modified for: Nextflow pipeline

Loads standardised embeddings and their true/predicted labels, reduces the
embedding space to 2D via t-SNE, and saves a two-panel scatter plot comparing
ground truth species labels against K-means cluster assignments.

Usage
-----
    python visualize_embeddings.py \\
        --embedding_file path/to/<base>_emb_stand.npy \\
        --label_file     path/to/<base>_labels.txt \\
        --clustering_file path/to/<base>_kmeans_results.csv \\
        --output_file    path/to/<base>_visualization.png

Arguments
---------
--embedding_file  : Standardised embedding matrix (.npy), shape (N_reads x embedding_dim).
                    Use _emb_stand.npy from calculate_embedding_for_tsv.py.
--label_file      : Labels file (.txt) from calculate_embedding_for_tsv.py.
                    First line must be the header "true_label".
--clustering_file : K-means results CSV from cluster_embeddings.py.
                    Must contain columns: true_label, predicted_labels.
--output_file     : Path to save the output PNG file (dpi=300).

Output
------
A two-panel PNG figure:
    Panel 1 — True Labels (Ground Truth): reads coloured by species.
    Panel 2 — K-means Clustering: reads coloured by predicted cluster ID.
Both panels share the same t-SNE 2D projection.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

def load_embeddings_from_npy(npy_file,label_file):
    """
    Loads standardised embedding matrix and true labels
    Parameters:
        npy_file: Path to .npy file containing the embedding matrix
        label_file: Path to labels .txt file
    """
    
    embeddings = np.load(npy_file)
    with open(label_file, 'r') as f:
        lines= f.read().strip().split('\n')
        labels= lines[1:]
    return embeddings, labels

def load_clustering_results(clustering_csv_file):
    """
    Loads true and predicted labels from a K-means results CSV
    Parameters:
        clustering_csv_file: Path to CSV produced by cluster_embeddings.py
    Returns:
        true_labels :True genome labels per read
        pred_labels : Predicted cluster IDs per read (best seed from K-means)d
    """
    # loads k-means clustering reuslts
    df = pd.read_csv(clustering_csv_file)
    true_labels = df['true_label'].values
    pred_labels = df['predicted_labels'].values

    return true_labels,pred_labels

# apply t-SNE to the embeddings to reduce their dimensionality to 2D
def dimreduct_TSNE(X=np.array):
    X_embedded = TSNE(n_components=2,learning_rate='auto',perplexity=100).fit_transform(X)
    return X_embedded

def create_visualization(embeddings_2d,true_labels,predicted_labels, output_file):
    '''
    embeddings_2d: 2D array after t-SNE
    true_labels: array with real species names
    predicted_labels: array with k-means clusterIDs
    output_file: String: "embeddings_visualization.png"
    '''
    unique_true = list(set(true_labels))
    # dict: name -> number 
    true_label2id = {l: i for i, l in enumerate(unique_true)}
    true_labels_numeric = np.array([true_label2id[l] for l in true_labels])
    # num of species
    num_clusters = len(unique_true)

    # create color maps 
    colors_true = plt.cm.viridis([i/(num_clusters-1) if num_clusters > 1 else 0 
                               for i in true_labels_numeric])
    colors_pred = plt.cm.viridis([i/(num_clusters-1) if num_clusters > 1 else 0 
                               for i in predicted_labels])

    fig, ax = plt.subplots(1, 2, figsize=(18, 6))
    
    # Plot 1: True labels 
    ax[0].scatter(embeddings_2d[:, 0], embeddings_2d[:, 1], c=colors_true, alpha=0.6)
    ax[0].set_title('True Labels (Ground Truth)', fontsize=14)
    ax[0].set_xlabel('t-SNE dimension 1')
    ax[0].set_ylabel('t-SNE dimension 2')
    
    # Plot 2: Predicted labels (K-means)
    ax[1].scatter(embeddings_2d[:, 0], embeddings_2d[:, 1], c=colors_pred, alpha=0.6)
    ax[1].set_title('K-means Clustering', fontsize=14)
    ax[1].set_xlabel('t-SNE dimension 1')
    ax[1].set_ylabel('t-SNE dimension 2')

    plt.tight_layout()
    
    # Save 
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Visualization saved to {output_file}")

def main(args):
    # Load standardized embeddings (for t-SNE)
    embeddings_std,labels = load_embeddings_from_npy(args.embedding_file,args.label_file)

    # Load clustering results
    true_labels, predicted_labels = load_clustering_results(args.clustering_file)
    
    # Perform t-SNE
    embeddings_2d = dimreduct_TSNE(embeddings_std)
    
    # Create visualization
    create_visualization(embeddings_2d, true_labels, predicted_labels, args.output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Visualize embeddings with t-sne')
    parser.add_argument('--embedding_file',type=str, required=True,help='Path to combined_embeddings_stand.npy(standardized)')
    parser.add_argument('--label_file',type=str, required=True,help='Path to combined_labels_stand.txt')
    parser.add_argument('--clustering_file',type=str, required=True,help='Path to kmeans_results.csv')
    parser.add_argument('--output_file',type=str, required=True,help='Path to save visualization PNG')
    args = parser.parse_args()
    main(args)
    
