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
    publishDir "${projectDir.parent}/data/species_embeddings",mode:'copy'

    input:
        file tsv_file
        path py_script
        path path_model_dir
    output:
        tuple   path("data/species_embeddings/${tsv_file.baseName}_emb.npy"),
                path("data/species_embeddings/${tsv_file.baseName}_labels.txt"),
                emit: normal
        tuple   path("data/species_embeddings/${tsv_file.baseName}_emb_stand.npy"),
                path("data/species_embeddings/${tsv_file.baseName}_labels.txt"),
                emit: standard
    script:
    """
    mkdir -p "data/species_embeddings" 

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
        path all_files
        path py_script
    output:
        tuple   path("combined_embeddings.npy"),
                path("combined_labels.txt"),
                emit: normal
        tuple   path("combined_embeddings_stand.npy"),
                path("combined_labels_stand.txt"),
                emit: standard
 
    script:
    """
    
    python3 ${py_script} \
            ${projectDir.parent}/data/species_embeddings \
            .
    """
}
process kmeans_clustering{
    publishDir "${projectDir.parent}/data/clustering_results",mode:'copy' 
    input:
        tuple path(combined_embeddings_npy), path(combined_labels_txt)
        path py_script
    output:
        path "kmeans_results.csv"
        path "kmeans_results_metrics.txt"
    script:
    """
    
    python3 ${py_script} \
            --embedding_file=${combined_embeddings_npy} \
            --label_file=${combined_labels_txt} \
            --output_file=kmeans_results.csv
    """
}

process visualization_embeddings{
    publishDir "${projectDir.parent}/data/visualizations",mode:'copy' 
    input:
        tuple path(combined_embeddings_stand_npy), path(combined_labels_stand_txt)
        path kmeans_results
        path py_script
    output:
        path "embeddings_visualization.png"
    script:
    """
    
    python3 ${py_script} \
            --embedding_file=${combined_embeddings_stand_npy} \
            --label_file=${combined_labels_stand_txt} \
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
    
    all_files = calculate_embeddings.out.normal
        .concat(calculate_embeddings.out.standard)
        .flatten()
        .collect()

    combine_embeddings(all_files, file("${projectDir}/combine_embeddings.py"))

    kmeans_clustering(combine_embeddings.out.normal, file("${projectDir}/cluster_embeddings.py"))

    visualization_embeddings(combine_embeddings.out.standard, kmeans_clustering.out[0], file("${projectDir}/visualize_embeddings.py"))
}