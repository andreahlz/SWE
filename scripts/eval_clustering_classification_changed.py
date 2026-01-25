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

from utils_changed import get_embedding

# wieder dekommentieren für ABgabe
#csv.field_size_limit(sys.maxsize)
# csv.field_size_limit(sys.maxsize)

# Set debug mode; true if args are hardcoded
DEBUG_MODE = False
device = "linux"
# muss ich noch am ende rauslöschen
if device == "windows":
    max_c_long = 2**31 - 1   # 2147483647
    csv.field_size_limit(max_c_long)

def dimreduct_TSNE(X=np.array):
    #todo: find the best value for perplexity, set perplexit to it. The larger the dataset, the higher the value should be 
    #run it several times, find the value, where kl divergence does not improve any more (KL divergence low indicates better results)
    #perplexity sets the effective number of neighbours that each point is attracted to
    #Kullback-Leibler (KL) divergence is a measure of the difference between two probability distributions (low and high dimensional)
    X_embedded = TSNE(n_components=2,learning_rate='auto',perplexity=2).fit_transform(X)
    return X_embedded

def main(args):
    model_list = args.model_list.split(",")
    for model in model_list:
        print(args.data_dir)
        for dataset in os.listdir(args.data_dir): #datasets. we only got one so far
            print(os.listdir(args.data_dir))
            embeddings_all = {}
            labels_data = {}
            colors_datapoints = {}
            for reads_mode in ["all_tsvs"]: 
                print(f"Start {model} {dataset} {reads_mode} clustering")
                data_file = os.path.join(args.data_dir, dataset, f"{reads_mode}.tsv")
                                
                with open(data_file, "r") as f:
                    reader = csv.reader(f, delimiter="\t")
                    data = list(reader)[1:]

                dna_sequences = [d[0] for d in data]
                labels = [d[1] for d in data]


                # convert labels to numeric values  
                label2id = {l: i for i, l in enumerate(set(labels))}
                labels = np.array([label2id[l] for l in labels])
                num_clusters = len(label2id)
                print(f"Get {len(dna_sequences)} sequences, {num_clusters} clusters")  
                #save true ground labels for the current reads mode in the dictionary
                #to colour the points in the visualization later
                labels_data[f"true_ground_{reads_mode}"] = labels

                # generate embedding
                embedding = get_embedding(dna_sequences, model, dataset, reads_mode, test_model_dir=args.test_model_dir)
                embedding_norm = normalize(embedding)
                #use embedding_standard for t-SNE
                embedding_standard = StandardScaler().fit_transform(embedding)
                #converts the values of the embedding into differences to the mean in standard deviations
                #reduce dimensionality of current embedding, save low dimensional embedding in dictionary
                embedding_reduced = dimreduct_TSNE(embedding_standard)
                embeddings_all[f"low_dim_{reads_mode}"] = embedding_reduced
                embeddings_all[f"embedding_{reads_mode}"] = embedding
                embeddings_all[f"embedding_standard_{reads_mode}"] = embedding_standard
                
                random_seeds = [0, 1, 2, 3, 4]
                #random_seeds = [0,1] #for quick testing, delete later
                kmeans_results = np.zeros([len(random_seeds), 4])
                lr_results = np.zeros([len(random_seeds), 5])
                
                for random_seed in random_seeds:
                    ### perform k-means clustering
                    kmeans = KMeans(n_clusters=num_clusters,
                                    random_state=random_seed,
                                    max_iter=1000,
                                    init="random",
                                    n_init=3)
                    kmeans.fit(embedding_norm)
                    preds_clustering = kmeans.labels_
                    #save predicted clustering labels for the current seed and reads mode in the dictionary
                    #todo compare the results for different random seeds later, take label with the most occurences
                    labels_data[f"kmeans_predicted_{reads_mode}"] = preds_clustering

                    purity = sklearn.metrics.homogeneity_score(labels, preds_clustering)
                    completeness = sklearn.metrics.completeness_score(labels, preds_clustering)
                    ari = sklearn.metrics.adjusted_rand_score(labels, preds_clustering)
                    nmi = sklearn.metrics.normalized_mutual_info_score(labels, preds_clustering)
                    
                    kmeans_results[random_seed] = np.array([purity, completeness, ari, nmi])
                    
                    print(f"Kmeans purity: {purity} completeness: {completeness} ari: {ari} nmi: {nmi}")
                    

                    # perform few-shot classification
                    results = []
                    
                    # generate different train/test splits for each random seed
                    np.random.seed(random_seed)
                    permutation = np.random.permutation(len(labels))
                    labels = labels[permutation]
                    embedding_norm = embedding_norm[permutation]
                    embedding_standard = embedding_standard[permutation]
                    
                    # for num_samples_per_class in [1, 2, 5, 10, 20]:
                    # #for num_samples_per_class in [2, 5, 20]: # for quick testing, delete later todo
                    #     is_train = np.zeros(len(labels))
                    #     is_test = np.zeros(len(labels))
                    #     for i in range(num_clusters):
                    #         idx = np.where(labels == i)[0]
                    #         is_train[idx[:num_samples_per_class]] = 1
                    #         is_test[idx[-80:]] = 1
                    #     is_train = is_train.astype(bool)
                    #     is_test = is_test.astype(bool)
                        
                    #     embedding_train = embedding_standard[is_train]
                    #     embedding_test = embedding_standard[is_test]

                    #     # 1. Logistic Regression
                    #     lr = LogisticRegression(random_state=random_seed, 
                    #                             max_iter=3000, 
                    #                             n_jobs=64,
                    #                             solver="lbfgs",
                    #                             penalty="l2",
                    #                             C=0.5)
                    #     lr.fit(embedding_train, labels[is_train])
                    #     preds_lr = lr.predict(embedding_test)
                    #     preds_train_lr = lr.predict(embedding_train)
                        
                    #     f1_train = sklearn.metrics.f1_score(labels[is_train], preds_train_lr, average="macro", zero_division=0)
                    #     loss_train = sklearn.metrics.log_loss(labels[is_train], lr.predict_proba(embedding_train))
                                        
                    #     f1 = sklearn.metrics.f1_score(labels[is_test], preds_lr, average="macro", zero_division=0)
                    #     recall = sklearn.metrics.recall_score(labels[is_test], preds_lr, average="macro", zero_division=0)
                    #     precision = sklearn.metrics.precision_score(labels[is_test], preds_lr, average="macro", zero_division=0)
                    #     accuracy = sklearn.metrics.accuracy_score(labels[is_test], preds_lr)
                    #     results.append(f1)
                    #     print(f"LR {num_samples_per_class}  train f1: {f1_train} loss: {loss_train} f1: {f1} recall: {recall} precision: {precision} accuracy: {accuracy}")
                    #todo: decomment later
                    #lr_results[random_seed] = np.array(results)
                
                kmeans_results = kmeans_results.mean(axis=0)
                print(f"Kmeans purity: {kmeans_results[0]} completeness: {kmeans_results[1]} ari: {kmeans_results[2]} nmi: {kmeans_results[3]}")
                
                lr_results = lr_results.mean(axis=0)
                lr_results = ", ".join([str(round(r, 6)) for r in lr_results])
                print(f"LR 1: {lr_results[0]} 2: {lr_results[1]} 5: {lr_results[2]} 10: {lr_results[3]} 20: {lr_results[4]}")
                print(lr_results)
                if f"{reads_mode}_labels_correct_predicted" not in colors_datapoints.keys():
                    print("test truth value calculation")
                    truth_value = [0 if labels_data[f"true_ground_{reads_mode}"][i] != labels_data[f"kmeans_predicted_{reads_mode}"][i] else 1 for i in range(len(labels_data[f"true_ground_{reads_mode}"]))]
                    colors_datapoints[f"{reads_mode}_labels_correct_predicted"] = ["green" if truth_value[i] == 1 else "red" for i in range(len(truth_value))]
            #visualization of embeddings
            #convert labels into colors
            colors = np.linspace(0,1,num_clusters)
            for l in labels_data.keys():
                labels_norm = [i/(num_clusters-1) for i in labels_data[l]]
                colors_datapoints[l] = plt.cm.viridis(labels_norm)
            print(embeddings_all.keys())
            fig,ax = plt.subplots(2,3, figsize=(18,12))
            for ind,reads_mode in enumerate(["all_tsvs"]):
                titles = ["Kmeans predicted clustering", "True ground clustering", "Kmeans predicted clustering with correct/incorrect labels"]
                datapoints_labels = [f"kmeans_predicted_{reads_mode}", f"true_ground_{reads_mode}", f"{reads_mode}_labels_correct_predicted"]
                for i in range(3):
                    dim1 = [embeddings_all[f"low_dim_{reads_mode}"][i][0] for i in range(len(embeddings_all[f"low_dim_{reads_mode}"]))]
                    dim2 = [embeddings_all[f"low_dim_{reads_mode}"][i][1] for i in range(len(embeddings_all[f"low_dim_{reads_mode}"]))]
                    ax[ind,i].scatter(dim1,dim2,c=colors_datapoints[datapoints_labels[i]])
                    if ind == 0:
                        ax[ind,i].set_title(titles[i])
            plt.savefig("embeddings_visualization.png")
            plt.show()
            #plot for both short and long reads
        print("test")

                    
                

if __name__ == "__main__":
    os.environ["TOKENIZERS_PARALLELISM"] = "false"
    if DEBUG_MODE:
        #muss ich noch ausbauen
        proj_root = Path(__file__).parent.parent
        args = argparse.Namespace(
            #todo: in bash skript download von DNABERTS starten, Pfad bei test_model_dir übergeben
            test_model_dir = "/home/barbara/evaluate/DNABERT-S",
            data_dir=f"{proj_root}/data/DNABERTS_input",
            model_list="test"         
        )
    elif DEBUG_MODE and device == "windows":
        args = argparse.Namespace(
            test_model_dir = Path(r"C:\Users\barba\Desktop\Bioinformatik\DNABERT-S\DNABERT-S"),
            data_dir=Path(r"C:\Users\barba\SWE\data\DNABERTS_input"),
            model_list="tnf, test",
            #added sample size
            sample_size=1
        )
    else:
        parser = argparse.ArgumentParser(description='Evaluate clustering')
        parser.add_argument('--test_model_dir', type=str, default="/root/trained_model", help='Directory to save trained models to test')
        parser.add_argument('--model_list', type=str, default="tnf, test", help='List of models to evaluate, separated by comma. Currently support [tnf, tnf-k, dnabert2, hyenadna, nt, test]')
        parser.add_argument('--data_dir', type=str, default="/root/data", help='Data directory')
        # added sample size
        parser.add_argument('--sample_number', type=int, default=1, help='Number of samples to evaluate')
        args = parser.parse_args()
    print(args)
    print(type(args))
    main(args) 
