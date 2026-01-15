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
    mkdir "data/csv_embeddings" -p

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
process combine_embeddings{
    publishDir "${projectDir.parent}/data/commbined_embedings",mode:'copy' 

    input:
        path csv_files
        path py_script
    output:
        path "data/commbined_embedings/combined_embeddings.csv"
        path "data/commbined_embedings/combined_embeddings_stand.csv"
    script:
    """
    mkdir "data/combined_embeddings" -p
    python3 ${py_script} ${projectDir.parent}/data/csv_embeddings .
        
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
}