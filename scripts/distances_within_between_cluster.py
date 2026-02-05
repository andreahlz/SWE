import pandas as pd
import argparse
import numpy as np
from scipy.spatial.distance import cdist

def calc_within_cluster_distances(centroid,cluster_data_points):
    distances = np.linalg.norm(cluster_data_points-centroid,axis=1)
    return distances,sum(distances)            

def main(args):
    labels = pd.read_csv(args.label_file)
    distances = []
    embedding = np.load(args.embedding_file) 
    data = pd.DataFrame(data=embedding,index=labels['true_label'])
    data=pd.DataFrame(embedding,index=labels['true_label'])
    labels_set = set(labels['true_label']) 
    mean_within_cluster_distances = []
    #calculate centroids
    centroids = {lab:data.loc[lab].mean() for lab in labels_set}
    for lab in labels_set:
        cur_distances,wcd = calc_within_cluster_distances(centroids[lab],data.loc[lab])
        distances.append(cur_distances)
        mean_within_cluster_distances.append(wcd/len(data.loc[lab]))
    #calculate between cluster distance using scipy (returns distance matrix)
    centroids_array= np.array(list(centroids.values()))
    between_dist = cdist(centroids_array,centroids_array)
    #save results as npz file
    np.savez(args.outfile,between = between_dist,mean_within = np.array(mean_within_cluster_distances),
             distances=np.array(distances).flatten(),true_label = labels['true_label'])

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Calculate distances within and between clusters')
    parser.add_argument('--embedding_file', type=str, help='Path to the embedding npy file')
    parser.add_argument('--label_file', type=str, help='Label file containing two columns. First for true label, second for predicted label')
    parser.add_argument('--outfile',type=str,default='calculated_distances.npz',help='Path tp npz file to save calculated distances.')
    args = parser.parse_args()
    main(args)