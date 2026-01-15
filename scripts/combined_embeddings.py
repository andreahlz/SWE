'''
combine embeddings: all individdual csv-embeddings of each species combined to 1
output of embedding calulations: 2 versions of embeddings (.csv and _stand.csv)
'''
import os
import pandas as pd
import sys


def combine_embeddings(csv_dir, output_dir):
    all_files = sorted(os.listdir(csv_dir))
    normal_paths = []
    stand_paths = []
    for file in all_files:
        file_path = os.path.join(csv_dir,file)

        #check if files with that directory excits
        if not os.path.isfile(file_path):
            continue

        #only csv files
        if not file.endswith(".csv"):
            continue

        if file.endswith("_stand.csv"):
            stand_paths.append(file_path)
        else:
            normal_paths.append(file_path)
    
    os.makedirs(output_dir,exist_ok=True) 
    dfs = []
    dfs_stand = []
    if normal_paths:
        for path in normal_paths:
            df = pd.read_csv(path,header=None)
            dfs.append(df)
        combined = pd.concat(dfs,ignore_index=True)
        combined.to_csv(os.path.join(output_dir, "combined_embeddings.csv"),header = False, index = False)
    
    if stand_paths:
        for path in stand_paths:
            df = pd.read_csv(path,header=None)
            dfs_stand.append(df)
        combined_stand = pd.concat(dfs_stand,ignore_index=True)
        combined_stand.to_csv(os.path.join(output_dir, "combined_embeddings_stand.csv"),header = False, index = False)


if __name__ == "__main__":
    if len(sys.arg) != 3:
        sys.stderr.write("Usage: python combine_embeddings.py <csv_dir> <output_dir>\n")
        sys.exit(2)
    combine_embeddings(sys.argv[1], sys.argv[2])