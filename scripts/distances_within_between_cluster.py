"""
distances_within_between_cluster.py

Calculates within- and between-cluster distances using Euclidean distance. 
Within-cluster distances are defined as the distance from each data point to its cluster centroid; 
between-cluster distances as centroid-to-centroid distances. 
Species labels are taken from the ground-truth clustering to ensure the same data points are used for both short and long read distance calculations. 
Results are saved to the .npz file specified by --outfile.

Usage
-----
    python distances_within_between_cluster.py \\
        --embedding_file path/to/<base>_emb.npy \\
        --label_file     path/to/<base>_labels.txt \\
        --outfile        path/to/<base>_cluster_distances.npz

Arguments
---------
--embedding_file: Path to raw embedding matrix (.npy), shape (N_reads x embedding_dim).
                  Use _emb.npy (not standardised) from calculate_embedding_for_tsv.py.
--label_file: Path to labels .txt file from calculate_embedding_for_tsv.py.
              First line must be the header "true_label" and is skipped.
--outfile: Path for the output .npz file. Default: calculated_distances.npz.
"""
import pandas as pd
import argparse
import numpy as np
from scipy.spatial.distance import cdist

def main(args):
    """
    Loads embeddings and labels, computes within- and between-cluster distances,
    and saves results to a .npz archive.

    Parameters
    args :
        Parsed arguments with attributes: embedding_file, label_file, outfile.
    """
    #get labels
    label_file = open(args.label_file,'r')
    labels = [i.rstrip("\n") for i in label_file.readlines()][1:]
    distances = []
    #save embedding and labels do dataframe
    embedding = np.load(args.embedding_file) 
    data=pd.DataFrame(embedding,index=labels)
    labels_set = set(labels) 
    #calculate centroids
    centroids = {lab:data.loc[[lab]].mean() for lab in labels_set}
    #initiate distance array, fill it with NaNs
    distances = np.empty(len(data))
    distances[:] = np.nan
    for lab,cent in centroids.items():
        #create mask, use mask to get data points
        pos = data.index == lab
        data_subset = data.loc[lab]
        #If there is more than one embedding for the current species
        if type(data_subset) == pd.DataFrame:
            distances[pos] = [np.linalg.norm(data_subset.values[i]-cent,axis=0) for i in range(len(data_subset))]
        #If the current species only contains one embedding
        else:
            distances[pos] = np.linalg.norm(data_subset-cent,axis = 0)
    #calculate between cluster distance using scipy (returns distance matrix)
    centroids_array= np.array(list(centroids.values()))
    between_dist = cdist(centroids_array,centroids_array)
    #save results as npz file
    np.savez(args.outfile,between = between_dist,
             distances=distances,true_label = labels)

if __name__=="__main__":
    # args = argparse.Namespace(
    #     embedding_file = '/home/barbara/SWE/data/distances_test/test_run_processed_emb.npy',
    #     label_file='/home/barbara/SWE/data/distances_test/test_run_processed_labels.txt',
    #     outfile = 'calculated_distances.npz'
    # )
    parser = argparse.ArgumentParser(description='Calculate distances within and between clusters')
    parser.add_argument('--embedding_file', type=str, help='Path to the embedding npy file')
    parser.add_argument('--label_file', type=str, help='Label file containing two columns. First for true label, second for predicted label')
    parser.add_argument('--outfile',type=str,default='calculated_distances.npz',help='Path tp npz file to save calculated distances.')
    args = parser.parse_args()
    main(args)