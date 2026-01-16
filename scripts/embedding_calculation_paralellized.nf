#!/usr/bin/env nextflow

process get_n_lines_of_tsv{
    publishDir "${projectDir.parent}/data/tsvs",mode:'copy'

    input:
        file tsv_file
    output:
        file "${tsv_file.baseName}_processed.tsv"
    script:
    """
    echo "${projectDir.parent}/data/tsvs"
    head ${tsv_file} -n 65000 > ${tsv_file.baseName}_processed.tsv
    """
}
process calculate_embeddings{
    publishDir "${projectDir.parent}/data/csv_embeddings",mode:'copy'

    input:
        file tsv_file
        path py_script
        path path_model_dir
    output:
        path "data/csv_embeddings/${tsv_file.baseName}.csv"
        path "data/csv_embeddings/${tsv_file.baseName}_stand.csv"
    script:
    """
    mkdir -p "data/csv_embeddings" 

    python3 ${py_script} \
            --tsv_file_path=${tsv_file} \
            --test_model_dir=${path_model_dir} \
            --model_list="test"
    """
}
/*
/missing: process to combine embeddings, TSNe on these,
/set argument path_model_dir such that it finds the model
/*/

process combine_embeddings{
    publishDir "${projectDir.parent}/data/combined_embeddings",mode:'copy' 

    input:
        path csv_files
        path py_script
    output:
        path "combined_embeddings.csv"
        path "combined_embeddings_stand.csv"
    script:
    """
    
    python3 ${py_script} \
            ${projectDir.parent}/data/csv_embeddings \
            .
    """
}
process kmeans_clustering{
    publishDir "${projectDir.parent}/data/clustering_results",mode:'copy' 
    input:
        path combined_embeddings
        path py_script
    output:
        path "kmeans_results.csv"
        path "kmeans_results_metrics.txt"
    script:
    """
    
    python3 ${py_script} \
            --embedding_file=${combined_embeddings} \
            --output_file=kmeans_results.csv
    """
}
process visualization_embeddings{
    publishDir "${projectDir.parent}/data/visualizations",mode:'copy' 
    input:
        path combined_embeddings_stand
        path kmeans_results
        path py_script
    output:
        path "embeddings_visualization.png"
    script:
    """
    
    python3 ${py_script} \
            --embedding_file=${combined_embeddings_stand} \
            --clustering_file=${kmeans_results} \
            --output_file=embeddings_visualization.png
    """
}




workflow{
    println projectDir
    data_dir = "${projectDir.parent}/data"
    println data_dir

    only_tsvs = file("${data_dir}/tsvs/*.tsv")
    println only_tsvs
    
    tsv_files = channel.fromPath("${data_dir}/tsvs/*.tsv")
    
    get_n_lines_of_tsv(tsv_files)
    
    calculate_embeddings(get_n_lines_of_tsv.out,"${projectDir}/calculate_embedding_for_tsv.py","/home/barbara/evaluate/DNABERT-S")
    
    all_csv_files = calculate_embeddings.out[0].collect()
    combine_embeddings(all_csv_files, file("${projectDir}/combine_embeddings.py"))

    kmeans_clustering(calculate_embeddings.out[0], file("${projectDir}/kmeans_clustering.py"))

    visualization_embeddings(combine_embeddings.out[1], kmeans_clustering.out[0], file("${projectDir}/visualize_embeddings.py"))
}