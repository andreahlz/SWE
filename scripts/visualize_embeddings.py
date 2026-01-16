"""
t-SNE visualization of embeddings with clustering results

Based on: eval_clustering_classification_changed.py 
Modified for: Nextflow pipeline
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

def load_embeddings_from_npy(npy_file,label_file):
    embeddings = np.load(npy_file)
    with open(label_file, 'r') as f:
        labels = f.read().strip().split('\n')
    return embeddings, labels

def load_clustering_results(clustering_csv_file):
    #loads k-means clustering reuslts
    df = pd.read_csv(clustering_csv_file)
    true_labels = df['true_label'].values
    pred_labels = df['predicted_labels'].values

    return true_labels,pred_labels

def dimreduct_TSNE(X=np.array):
    #todo: find the best value for perplexity, set perplexit to it. The larger the dataset, the higher the value should be 
    #run it several times, find the value, where kl divergence does not improve any more (KL divergence low indicates better results)
    #perplexity sets the effective number of neighbours that each point is attracted to
    #Kullback-Leibler (KL) divergence is a measure of the difference between two probability distributions (low and high dimensional)
    X_embedded = TSNE(n_components=2,learning_rate='auto',perplexity=30).fit_transform(X)
    return X_embedded

def create_visualization(embeddings_2d,true_labels,predicted_labels, output_file):
    '''
    embeddings_2d: 2D array after t-SNE
    true_labels: array with real species names
    predicted_labels: array with k-means clusterIDs
    output_file: String: "embeddings_visualization.png"
    '''
    unique_true = list(set(true_labels))
    # dic: name -> number 
    true_label2id = {l: i for i, l in enumerate(unique_true)}
    true_labels_numeric = np.array([true_label2id[l] for l in true_labels])
    #num of species
    num_clusters = len(unique_true)

    #creates color maps 
    colors_true = plt.cm.viridis([i/(num_clusters-1) if num_clusters > 1 else 0 
                               for i in true_labels_numeric])
    colors_pred = plt.cm.viridis([i/(num_clusters-1) if num_clusters > 1 else 0 
                               for i in predicted_labels])
    
    correct = true_labels_numeric == predicted_labels
    colors_correct = ['green' if c else 'red' for c in correct]

    fig, ax = plt.subplots(1, 3, figsize=(18, 6))
    
    # Plot 1: True labels 
    ax[0].scatter(embeddings_2d[:, 0], embeddings_2d[:, 1], c=colors_true, alpha=0.6)
    ax[0].set_title('True Labels (Ground Truth)', fontsize=14)
    ax[0].set_xlabel('t-SNE dimension 1')
    ax[0].set_ylabel('t-SNE dimension 2')
    
    # Plot 2: Predicted labels (K-means)
    ax[1].scatter(embeddings_2d[:, 0], embeddings_2d[:, 1], c=colors_pred, alpha=0.6)
    ax[1].set_title('K-means Predicted Clustering', fontsize=14)
    ax[1].set_xlabel('t-SNE dimension 1')
    ax[1].set_ylabel('t-SNE dimension 2')
    
    # Plot 3: Correct/Incorrect
    ax[2].scatter(embeddings_2d[:, 0], embeddings_2d[:, 1], c=colors_correct, alpha=0.6)
    ax[2].set_title('Clustering Quality (Green=Correct, Red=Incorrect)', fontsize=14)
    ax[2].set_xlabel('t-SNE dimension 1')
    ax[2].set_ylabel('t-SNE dimension 2')
    
    plt.tight_layout()
    
    # Save 
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Visualization saved to {output_file}")
    
    # Calculate accuracy for info
    accuracy = np.sum(correct) / len(correct)

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
    