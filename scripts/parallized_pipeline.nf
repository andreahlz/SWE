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
params.dataset_name = "mock_test"
params.data_dir_name = "data_test"
params.input_fastq = "${projectDir.parent}/data/simulated_reads/short_reads"
params.output_base = "${projectDir.parent}/${params.data_dir_name}"
params.config = "${projectDir.parent}/config.env"
params.input_txt= "${projectDir.parent}/${params.data_dir_name}/input/input.txt"


//python scripts

params.py_fastq_to_tsv = "${projectDir}/fastq_to_tsv.py"
params.py_embedding_calc = "${projectDir}/calculate_embedding_for_tsv_02.py"
params.py_kmeans = "${projectDir}/cluster_embeddings.py"
params.py_visualization = "${projectDir}/visualize_embeddings.py"
params.py_download = "${projectDir}/download_genome.py"
params.py_process_fna = "${projectDir}/process_genomes.py"

//simulation parameters
params.mean_length=1000
/*
params.mean_identity=95
params.identity_stdev=5

*/ 
params.length_stdev=500
params.indentity="nanopore2023"
params.seed=42

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
    publishDir "${params.output_base}/short_reads/simulated_reads", mode: 'copy'
    input:
        path fasta_file
        path input_txt
    output:
        tuple val('short_reads'), path("${fasta_file.baseName}.fastq")
        
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
process long_reads {
    publishDir "${params.output_base}/long_reads/simulated_reads", mode: 'copy'
    
    input:
        path fasta_file
        path input_txt
    
    output:
        tuple val('long_reads'), path("${fasta_file.baseName}.fastq")
        
    
    shell:
    '''
    #!/bin/bash
    set -euo pipefail
    
    output_dir="."
    file_name="!{fasta_file.baseName}"
    temp_dir="${output_dir}/temp"
    mkdir -p "${temp_dir}"

    # Read abundances
    declare -A abundances
    while IFS=$'\\t' read -r species abundance; do
        [[ "$species" == "species" ]] && continue
        abundances["$species"]=$abundance
    done < "!{input_txt}"

    # Split multifasta
    echo "Splitting multifasta into individual genomes..."
    awk -v outdir="${temp_dir}" '
    /^>/ {
        if (seq != "" && filename != "") {
            print seq > filename
            close(filename)
        }
        seq = ""
        split($0, a, " ")
        gsub(/^>/, "", a[1])
        filename = outdir "/" a[1] ".fna"
        print $0 > filename
        next
    }
    {
        seq = seq $0
    }
    END {
        if (seq != "" && filename != "") {
            print seq > filename
            close(filename)
        }
    }' "!{fasta_file}"

    echo "Created genome files:"
    ls -lh "${temp_dir}/"*.fna

    # Simulate reads
    echo ""
    echo "Simulating long reads with Badread..."
    for genome_file in "${temp_dir}"/*.fna; do
        [[ ! -f "$genome_file" ]] && continue
        
        genome_name=$(basename "$genome_file" .fna)
        abundance=${abundances[$genome_name]:-0.01}
        coverage=$(echo "$abundance * 10" | bc -l | awk '{printf "%.2f", $0}')
        
        echo "  ${genome_name}: abundance=${abundance}, coverage=${coverage}x"
        
        badread simulate --reference "${genome_file}" --quantity "${coverage}x" --length !{params.mean_length},!{params.length_stdev} --seed !{params.seed} > "${temp_dir}/${genome_name}.fastq"
    done
    
    echo ""
    echo "Created FASTQ files:"
    ls -lh "${temp_dir}"/*.fastq
    
    echo ""
    echo "Combining reads..."
    cat "${temp_dir}"/*.fastq > "${output_dir}/${file_name}.fastq"

    echo "Cleaning up..."
    rm -rf "${temp_dir}"
    
    echo ""
    echo "✓ Done! Output: ${file_name}.fastq"
    '''
}
//Convert FASTQ to TSV format
process converting_reads_format {
  publishDir "${params.output_base}/${read_type}/simulated_reads", mode: 'copy'

  input:
    tuple val(read_type), path(fastq_file)

  output:
    tuple val(read_type), path("${fastq_file.simpleName}.tsv")

  script:
  """
  python3 "${params.py_fastq_to_tsv}" \
    --input "${fastq_file}" \
    --output "${fastq_file.simpleName}.tsv"
  """
}
// only for testing 
process get_n_lines_of_tsv {
    publishDir "${params.output_base}/${read_type}/tsv_processed", mode: 'copy'
    
    input:
        tuple val(read_type),path(tsv_file)
    
    output:
        tuple val(read_type),path("${tsv_file.baseName}_processed.tsv")
    
    script:
    """
    head -n 10 ${tsv_file} > ${tsv_file.baseName}_processed.tsv
    echo "Reduced from \$(wc -l < ${tsv_file}) to 10 lines"
    """
}
//Calculate embeddings
process calculate_embeddings {
  publishDir "${params.output_base}/${read_type}/embeddings", mode: 'copy'

  input:
    tuple val(read_type), path(tsv_file)

  output:
    tuple val(read_type), path("${tsv_file.simpleName}_emb.npy")
    tuple val(read_type), path("${tsv_file.simpleName}_emb_stand.npy")
    tuple val(read_type), path("${tsv_file.simpleName}_labels.txt")

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
  publishDir "${params.output_base}/${read_type}/clustering", mode: 'copy'

  input:
    tuple val(read_type), path(embedding_npy), path(labels_txt)

  output:
    tuple val(read_type), path("${embedding_npy.simpleName}_kmeans_results.csv")
    tuple val(read_type), path("${embedding_npy.simpleName}_kmeans_results_metrics.txt")

  script:
  """
  python3 "${params.py_kmeans}" \
    --embedding_file "${embedding_npy}" \
    --label_file "${labels_txt}" \
    --output_file "${embedding_npy.simpleName}_kmeans_results.csv"
    """
}
process visualization_embeddings{
    publishDir "${params.output_base}/${read_type}/visualizations",mode:'copy' 
    input:
        tuple val(read_type), path(embeddings_stand_npy), path(labels_txt), path(kmeans_results)
    output:
        tuple val(read_type), path("${embeddings_stand_npy.simpleName}.png")
    script:
    """
    
    python3 ${params.py_visualization} \
            --embedding_file ${embeddings_stand_npy} \
            --label_file ${labels_txt} \
            --clustering_file ${kmeans_results} \
            --output_file "${embeddings_stand_npy.simpleName}.png"
    """
}
workflow process_read_pipeline {
    take: //declares the inputs of a named workflow
        fastq_tuple // tuple (read_type, fastqfile)
    main: 
    //FASTQ to TSV
    tsv_ch = converting_reads_format(fastq_tuple)

    //TSV shorter (for testing)
    tsv_processed = get_n_lines_of_tsv(tsv_ch)

    //Calculate embeddings
    emb_outputs = calculate_embeddings(tsv_processed)
    // emb_outputs[0] = tuple(read_type, emb.npy)
    // emb_outputs[1] = tuple(read_type, emb_stand.npy)
    // emb_outputs[2] = tuple(read_type, labels.txt)

    // keans clustering
    kmenas_input = emb_outputs[0].join(emb_outputs[2])
    //kmenas_input=    
    // tuple value(read_type), path(embedding_npy)
    // tuple value(read_type), path(labels_txt)
    
    kmeans_out=kmeans_clustering(kmenas_input)
    // tuple value(read_type), path("${embedding_npy.simpleName}_kmeans_results.csv")
    // tuple value(read_type), path("${embedding_npy.simpleName}_kmeans_results_metrics.txt")

    // Visualization
    //tuple val(read_type), path(embeddings_stand_npy), path(labels_txt), path(kmeans_results)
    viz_input=emb_outputs[1].join(emb_outputs[2]).join(kmeans_out[0])
    visualization_embeddings(viz_input)
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

    short_ch = short_reads(processed_fastas[1], input_ch)
    long_ch  = long_reads(processed_fastas[1], input_ch)

    reads_ch = short_ch.mix(long_ch)
    process_read_pipeline(reads_ch)
}
