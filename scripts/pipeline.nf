#!/usr/bin/env nextflow

/*
 * Pipeline: DNABERT-S Embedding Analysis Pipeline
 * Description: Converts FASTQ reads to TSV, calculates embeddings, performs clustering, and visualizes results
 * 
 * Workflow:
 * 1. converting_reads_format: FASTQ → TSV conversion
 * 2. calculate_embeddings: Generate DNABERT-S embeddings
 * 3. kmeans_clustering: Cluster embeddings using K-means
 * 4. visualization: Create t-SNE visualization of clusters
 */ 

// Parameters -> will be written in config
params.read_type = "short_reads"
params.input_fastq = "${projectDir.parent}/data/simulated_reads/short_reads"
params.model_list = "short_reads"
params.output_base = "${projectDir.parent}/data_test/${params.read_type}"

//python scripts
params.py_fastq_to_tsv = "${projectDir}/fastq_to_tsv.py"
params.py_embedding_calc = "${projectDir}/calculate_embedding_for_tsv_02.py"
params.py_kmeans = "${projectDir}/cluster_embeddings.py"
params.py_visualization = "${projectDir}/visualize_embeddings.py"


// DNABERT-S local directory (used with model_list "test")
params.test_model_dir = "${projectDir.parent}/DNABERT-S"

params.data_dir = "."

//combining two fastq outputfiles
process combine_read_output {
    publishDir "${params.output_base}/simulated_reads", mode: 'copy'
    
    input:
        tuple val(base_name), path(reads)
    
    output:
       path "${base_name}.fastq"
       
    script:

    """
    cat ${reads} > ${base_name}.fastq
    """
}

//Convert FASTQ to TSV format
process converting_reads_format {
  publishDir "${params.output_base}/simulated_reads", mode: 'copy'

  input:
    path fastq_file

  output:
    path "${fastq_file.simpleName}.tsv"

  script:
  """
  python3 "${params.py_fastq_to_tsv}" \
    --input "${fastq_file}" \
    --output "${fastq_file.simpleName}.tsv"
  """
}
//Calculate embeddings
process calculate_embeddings {
  publishDir "${params.output_base}/embeddings", mode: 'copy'

  input:
    path tsv_file

  output:
    path "${tsv_file.simpleName}_emb.npy"
    path "${tsv_file.simpleName}_emb_stand.npy"
    path "${tsv_file.simpleName}_labels.txt"

  script:
  """
  python3 "${params.py_embedding_calc}" \
    --tsv_file_path="${tsv_file}" \
    --model_list="test" \
    --test_model_dir="${params.test_model_dir}" \
    --data_dir="${params.data_dir}" \
    --out_dir .
    """
}
process checkout_all_outputs {
    echo true
    input:
        path emb_file
        path stand_file
        path labels_file

    script:
    """
    echo "=== Checking embeddings ==="
    python3 "${projectDir}/checkout_npy.py" "${emb_file}"
    
    echo ""
    echo "=== Checking standardized embeddings ==="
    python3 "${projectDir}/checkout_npy.py" "${stand_file}"
    
    echo ""
    echo "=== Checking labels ==="
    wc -l "${labels_file}"
    head -n 5 "${labels_file}"
    """
}
//Calculate kmean clusters
process kmeans_clustering {
  publishDir "${params.output_base}/clustering", mode: 'copy'

  input:
    path embedding_npy
    path labels_txt

  output:
    path "${embedding_npy.simpleName}_kmeans_results.csv"
    path "${embedding_npy.simpleName}_kmeans_results_metrics.txt"

  script:
  """
  python3 "${params.py_kmeans}" \
    --embedding_file "${embedding_npy}" \
    --label_file "${labels_txt}" \
    --output_file "${embedding_npy.simpleName}_kmeans_results.csv"
    """
}
process visualization_embeddings{
    publishDir "${params.output_base}/visualizations",mode:'copy' 
    input:
        path embeddings_stand_npy
        path labels_txt
        path kmeans_results
    output:
        path "${embeddings_stand_npy.simpleName}.png"
    script:
    """
    
    python3 ${params.py_visualization} \
            --embedding_file ${embeddings_stand_npy} \
            --label_file ${labels_txt} \
            --clustering_file ${kmeans_results} \
            --output_file "${embeddings_stand_npy.simpleName}.png"
    """
}

workflow{
    println "current directory ${projectDir}"
    println "Input FASTQ: ${params.input_fastq}"
    println "Output dir:  ${params.output_base}"
    println "FASTQ->TSV script: ${params.py_fastq_to_tsv}"
    println "Embedding script:  ${params.py_embedding_calc}"
    println "DNABERT-S dir:     ${params.test_model_dir}"
    println projectDir
    data_dir = "${projectDir.parent}/data"
    reads = Channel.fromFilePairs("${params.input_fastq}/*_R{1,2}.fastq")
    
    combined = combine_read_output(reads)
    tsv_ch = converting_reads_format(combined)
    emb_outputs=calculate_embeddings(tsv_ch)
    checkout_all_outputs(emb_outputs[0], emb_outputs[1], emb_outputs[2])
    kmeans_out=kmeans_clustering(emb_outputs[0],emb_outputs[2])
    visualization_embeddings(emb_outputs[1],emb_outputs[2],kmeans_out[0])
}
