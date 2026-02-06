import pandas as pd
import argparse
import numpy as np
from scipy.spatial.distance import cdist

def calc_within_cluster_distances(centroid,cluster_data_points):
    if type(cluster_data_points) == pd.DataFrame:
        distances = np.linalg.norm(cluster_data_points-centroid,axis=1)
        return distances,sum(distances)  
    else:
        distances = np.array([np.linalg.norm(cluster_data_points - centroid)])
        return distances,distances.sum()

def main(args):
    labels = pd.read_csv(args.label_file)
    label_file = open(args.label_file,'r')
    labels = [i.rstrip("\n") for i in label_file.readlines()][1:]
    distances = []
    embedding = np.load(args.embedding_file) 
    embedding_transposed = embedding.transpose()
    data=pd.DataFrame(embedding,index=labels)
    labels_set = set(labels) 
    mean_within_cluster_distances = []
    #calculate centroids
    centroids = {lab:data.loc[[lab]].mean() for lab in labels_set}
    distances = np.empty(len(data))
    distances[:] = np.nan
    for lab,cent in centroids.items():
        pos = data.index == lab
        data_subset = data.loc[lab]
        if type(data_subset) == pd.DataFrame:
            distances[pos] = [np.linalg.norm(data_subset.values[i]-cent,axis=0) for i in range(len(data_subset))]
        else:
            distances[pos] = np.linalg.norm(data_subset-cent,axis = 0)
    #calculate between cluster distance using scipy (returns distance matrix)
    centroids_array= np.array(list(centroids.values()))
    between_dist = cdist(centroids_array,centroids_array)
    #save results as npz file
    np.savez(args.outfile,between = between_dist,
             distances=distances,true_label = labels)

if __name__=="__main__":
    args = argparse.Namespace(
        embedding_file = '/home/barbara/SWE/data/distances_test/test_run_processed_emb.npy',
        label_file='/home/barbara/SWE/data/distances_test/test_run_processed_labels.txt',
        outfile = 'calculated_distances.npz'
    )
    # parser = argparse.ArgumentParser(description='Calculate distances within and between clusters')
    # parser.add_argument('--embedding_file', type=str, help='Path to the embedding npy file')
    # parser.add_argument('--label_file', type=str, help='Label file containing two columns. First for true label, second for predicted label')
    # parser.add_argument('--outfile',type=str,default='calculated_distances.npz',help='Path tp npz file to save calculated distances.')
    # args = parser.parse_args()
    main(args)