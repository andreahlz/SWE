#!/usr/bin/env nextflow

process calculate_embeddings{
    publishDir "data/embeddings/${tsv_file.baseName}",mode:'copy'

    input:
        file tsv_file
        path py_script
        path path_model_dir
    output:
        path "data/csv_embeddings/${tsv_file.baseName}.csv"
        path "data/csv_embeddings/${tsv_file.baseName}_stand.csv"
    script:
    """
    echo "${tsv_file}"
    mkdir "data/embeddings/${tsv_file.baseName}" -p
    mkdir "data/csv_embeddings"
    python3 ${py_script}\
            --tsv_file_path=${tsv_file}\
            --test_model_dir=${path_model_dir}\
            --model_list="test"
    """
}
/*
/missing: process to combine embeddings, TSNe on these,
/set argument path_model_dir such that it finds the model
/*/
workflow{
    println projectDir
    data_dir = "${projectDir.parent}/data"
    println data_dir
    only_tsvs = file("${data_dir}/tsvs/*.tsv")
    println only_tsvs
    tsv_files = channel.fromPath("${data_dir}/tsvs/*.tsv")
    calculate_embeddings(tsv_files,"${projectDir}/calculate_embedding_for_tsv.py","/home/barbara/evaluate/DNABERT-S")

}