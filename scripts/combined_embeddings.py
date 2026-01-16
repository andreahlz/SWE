'''
combine embeddings: all individdual csv-embeddings of each species combined to 1
output of embedding calulations: 2 versions of embeddings (.csv and _stand.csv)
'''
import os
import numpy as np
import sys


def combine_embeddings(npy_dir, output_dir):
    all_files = sorted(os.listdir(npy_dir))
    normal_files = []
    stand_files = []
    label_files = []

    for file in all_files:
        file_path = os.path.join(npy_dir,file)

        #check if files with that directory excits
        if not os.path.isfile(file_path):
            continue
        

        #only csv files
        if file.endswith("_emb_stand.npy"):
            stand_files.append(file_path)
        elif file.endswith("_emb.npy"):
            normal_files.append(file_path)

        elif file.endswith("_labels.txt"):
            label_files.append(file_path)
    
    os.makedirs(output_dir,exist_ok=True) 


    if normal_files:
        embeddings = []
        labels = []
        for emb_file in normal_files:
            emb = np.load(emb_file)
            embeddings.append(emb)

            base_name = os.path.basename(emb_file).replace('_emb.npy', '')
            label_file = os.path.join(npy_dir, f"{base_name}_labels.txt")
            with open(label_file, 'r') as f:
                file_labels = f.read().strip().split('\n')
            labels.extend(file_labels)

        combined_emb = np.vstack(embeddings)
        np.save(os.path.join(output_dir, "combined_embeddings.npy"), combined_emb)
        with open(os.path.join(output_dir, "combined_labels.txt"), 'w') as f:
            f.write('\n'.join(labels))

    if stand_files:
        embeddings_stand = []
        labels_stand = []
        for path in stand_files:
            emb = np.load(path)
            embeddings_stand.append(emb)

            base_name = os.path.basename(path).replace('_emb_stand.npy', '')
            label_file = os.path.join(npy_dir, f"{base_name}_labels.txt")
            if os.path.exists(label_file):
                with open(label_file, 'r') as f:
                    file_labels = f.read().strip().split('\n')
                labels_stand.extend(file_labels)
            else:
                print(f"  WARNING: No labels found for {path}")
                labels_stand.extend(['unknown'] * len(emb))
            
            print(f"  Added {emb.shape[0]} sequences from {base_name}")
        
        combined_emb_stand = np.vstack(embeddings_stand)
        
        output_stand_file = os.path.join(output_dir, "combined_embeddings_stand.npy")
        np.save(output_stand_file, combined_emb_stand)  

        output_label_stand_file = os.path.join(output_dir, "combined_labels_stand.txt")
        with open(output_label_stand_file, 'w') as f:
            f.write('\n'.join(labels_stand))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: python combine_embeddings.py <npy_dir> <output_dir>\n")
        sys.exit(2)
    combine_embeddings(sys.argv[1], sys.argv[2])