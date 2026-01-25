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

 usage: nextflow run pipeline.nf
 */ 

// Parameters -> will be written in config
params.read_type = "short_reads"
params.dataset_name = "mock_test"
params.data_dir_name = "data"
params.input_fastq = "${projectDir.parent}/data/simulated_reads/short_reads"
params.output_base = "${projectDir.parent}/${params.data_dir_name}/${params.read_type}"
params.config = "${projectDir.parent}/config.env"
params.input_txt= "${projectDir.parent}/${params.data_dir_name}/input/input.txt"


//python scripts

params.py_fastq_to_tsv = "${projectDir}/fastq_to_tsv.py"
params.py_embedding_calc = "${projectDir}/calculate_embedding_for_tsv_02.py"
params.py_kmeans = "${projectDir}/cluster_embeddings.py"
params.py_visualization = "${projectDir}/visualize_embeddings.py"
params.py_download = "${projectDir}/download_genome.py"
params.py_process_fna = "${projectDir}/process_genomes.py"

// DNABERT-S local directory (used with model_list "test")
params.test_model_dir = "${projectDir.parent}/DNABERT-S"

params.data_dir = "."

process download_fasta{
    tag "ncbi-genomes"
    publishDir "${params.output_base}/genomes/fasta_files", mode: 'copy'
    input:
        path input_txt
        
    output:
       path "*.fna"
       
    script:
    """
        python3 ${params.py_download} \
        --input_file "${input_txt}" \
        --output_dir . \
        --timeout 300
"""
}

process process_fasta{
    publishDir "${params.output_base}/genomes", mode: 'copy'
    input:
        path fna_files
        path input_txt
        
    output:
        path "processed/*fna"
        path "${params.dataset_name}.fna" 
        path "header_info.csv"
    script:
    """
    mkdir -p raw processed
    mv *.fna raw/

    python3 ${params.py_process_fna} \
        --input_file "${input_txt}" \
        --single_fastas_dir raw \
        --output_dir processed  \
        --csv_output "header_info.csv" \
        --multifasta_output "${params.dataset_name}.fna"
"""
}
//
/*
* Apply Insilicoseq to fasta file
*/ 
process short_reads {
    publishDir "${params.output_base}/simulated_reads", mode: 'copy'
    input:
        path fasta_file
        path input_txt
    output:
        path "${fasta_file.baseName}.fastq"
    script:
        """
        tail -n +2 "${input_txt}" > abundance.txt
        echo 'starting read generation for ${fasta_file.baseName}'
        echo 'fasta_file name = ${fasta_file}'
        iss generate \
            --genomes ${fasta_file} \
            --abundance_file abundance.txt \
            --model MiSeq \
            --n_reads 20 \
            --output temp_reads
        
        cat temp_reads_R1.fastq temp_reads_R2.fastq > ${fasta_file.baseName}.fastq
        """
}
/*
process long_reads {
    publishDir "${params.output_base}/simulated_reads", mode: 'copy'
    input:
        path fasta_file //multifasta
        path input_txt
    output:
        path "${fasta_file.baseName}.fastq"
    script:
        """
        mkdir -p temp_genomes

        echo "Splitting multifasta into individual genomes..."
        awk '/^>/ {
            if (seq) print seq; 
            seq=""; 
            split($0,a," "); 
            gsub(/^>/,"",a[1]); 
            filename=a[1]".fna"; 
            print $0 > "'${temp_genomes}'/"filename; 
            next
        } 
        {seq=seq$0} 
        END {if (seq) print seq}' "${fasta_file}"


        tail -n +2 "${input_txt}" > abundance.txt

        echo 'starting read generation for ${fasta_file.baseName}'
        echo 'fasta_file name = ${fasta_file}'
        iss generate \
            --genomes ${fasta_file} \
            --abundance_file abundance.txt \
            --model MiSeq \
            --n_reads 20 \
            --output temp_reads
        
        cat temp_reads_R1.fastq temp_reads_R2.fastq > ${fasta_file.baseName}.fastq
        """
}
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
*/ 
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

    input_ch = Channel.fromPath(params.input_txt)
    fasta_files = download_fasta(input_ch)
    processed_fastas = process_fasta(fasta_files.collect(), input_ch)
    combined_reads = short_reads(processed_fastas[1], input_ch)
    //reads = Channel.fromFilePairs("${params.input_fastq}/*_R{1,2}.fastq")
    
    tsv_ch = converting_reads_format(combined_reads)
    emb_outputs=calculate_embeddings(tsv_ch)
    checkout_all_outputs(emb_outputs[0], emb_outputs[1], emb_outputs[2])
    kmeans_out=kmeans_clustering(emb_outputs[0],emb_outputs[2])
    visualization_embeddings(emb_outputs[1],emb_outputs[2],kmeans_out[0])
}
