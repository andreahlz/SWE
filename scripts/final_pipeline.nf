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

process calculate_coverage_for_short_reads {
    publishDir "${params.output_base}/coverage", mode: 'copy'
    
    input:
        path processed_fna
        path input_txt
        
    output:
        path "coverage_short.txt"
    
    script:
    """
    mkdir -p processed
    mv *.fna processed/
    
    python3 ${params.py_coverage_calc} \
        --input_file "${input_txt}" \
        --newheader_fastas_dir processed \
        --target_coverage "${params.target_coverage}" \
        --short_read_length "${params.short_read_length}" \
        --txt_short_output coverage_short.txt 
    """
}

process short_reads {
    publishDir "${params.output_base}/short_reads/simulated_reads", mode: 'copy'
    input:
        path fasta_file
        path input_txt
        path coverage_short
    output:
        tuple val('short_reads'), path("${fasta_file.baseName}.fastq")
        
    script:
        """
        n_reads=\$(cat ${coverage_short})
        tail -n +2 "${input_txt}" > abundance.txt
        echo 'starting read generation for ${fasta_file.baseName}'
        echo 'fasta_file name = ${fasta_file}'
        iss generate \
            --genomes ${fasta_file} \
            --abundance_file abundance.txt \
            --model ${params.model_short} \
            --n_reads \${n_reads} \
            --seed ${params.seeds_sim} \
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
    set -euo pipefail
    
    output_dir="."
    file_name="!{fasta_file.baseName}"
    temp_dir="${output_dir}/temp"
    mkdir -p "${temp_dir}"

    # Read abundances
    #declare -A abundances
    #while IFS=$'\\t' read -r species abundance; do
    #   [[ "$species" == "species" ]] && continue
    #    abundances["$species"]=$abundance
    #done < "!{input_txt}"
    

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
        
        #abundance=${abundances[$genome_name]:-0.01}
        #n_reads=$(grep "^${genome_name}" ${coverage_long} | awk '{gsub(/,/, "", $2); print $2}')
        
        
        badread simulate \
            --reference "${genome_file}" \
            --quantity "!{params.target_coverage}x" \
            --length !{params.mean_length_long},!{params.length_stdev_long} \
            --identity !{params.identity} \
            --error_model !{params.error_model} \
            --seed !{params.seeds_sim} \
            --junk_reads 0 \
            --random_reads 0 \
            --chimeras 0 \
            > "${temp_dir}/${genome_name}.fastq"
    done
    
    echo ""
    echo "Created FASTQ files:"
    ls -lh "${temp_dir}"/*.fastq
    
    #echo ""
    #echo "Combining reads..."
    cat "${temp_dir}"/*.fastq > "${output_dir}/${file_name}.fastq"
    
    

    echo "Cleaning up..."
    rm -rf "${temp_dir}"
    
    echo ""
    echo "Done! Output: ${file_name}.fastq"
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
    head -n ${params.n_lines_test} ${tsv_file} > ${tsv_file.baseName}_processed.tsv
    echo "Reduced from \$(wc -l < ${tsv_file}) to ${params.n_lines_test} lines"
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
    debug true
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
process cluster_distance_calculations{
    publishDir "${params.output_base}/${read_type}/clustering",mode:'copy' 
    input:
        tuple val(read_type), path(embeddings_npy), path(labels_txt)
    output:
        tuple val(read_type), path("${embeddings_npy.simpleName}_${read_type}_cluster_distances.npz")
    script:
    """
    
    python3 ${params.py_cluster_dist} \
            --embedding_file ${embeddings_npy} \
            --label_file ${labels_txt} \
            --outfile "${embeddings_npy.simpleName}_${read_type}_cluster_distances.npz"
    """
}
process cluster_distance_visualization{
    publishDir "${params.output_base}/comparison",mode:'copy' 
    input:
        path short_reads_distances_npz
        path long_reads_distances_npz
    output:
        path "cluster_comparison.png"
    script:
    """
    
    python3 ${params.py_visualization_cluster_dist} \
            --file_1 ${short_reads_distances_npz} \
            --file_2 ${long_reads_distances_npz} \
            --result_png "cluster_comparison.png"
    """
}
workflow process_read_pipeline {
    take: //declares the inputs of a named workflow
        fastq_tuple // tuple (read_type, fastqfile)
    main: 
    //FASTQ to TSV
    tsv_ch = converting_reads_format(fastq_tuple)

    tsv_for_embddings = tsv_ch
    //TSV shorter (for testing)
    if( params.use_get_n_lines ) {
        tsv_for_embddings = get_n_lines_of_tsv(tsv_ch)
    }

    //Calculate embeddings
    emb_outputs = calculate_embeddings(tsv_for_embddings)
    // emb_outputs[0] = tuple(read_type, emb.npy)
    // emb_outputs[1] = tuple(read_type, emb_stand.npy)
    // emb_outputs[2] = tuple(read_type, labels.txt)

    // keans clustering
    kmeans_input = emb_outputs[0].join(emb_outputs[2])
    //kmeans_input=    
    // tuple value(read_type), path(embedding_npy)
    // tuple value(read_type), path(labels_txt)
    
    kmeans_out=kmeans_clustering(kmeans_input)
    // tuple value(read_type), path("${embedding_npy.simpleName}_kmeans_results.csv")
    // tuple value(read_type), path("${embedding_npy.simpleName}_kmeans_results_metrics.txt")

    // Visualization
    //tuple val(read_type), path(embeddings_stand_npy), path(labels_txt), path(kmeans_results)
    viz_input=emb_outputs[1].join(emb_outputs[2]).join(kmeans_out[0])
    visualization_embeddings(viz_input)
    calc_dist_cluster_input = emb_outputs[0].join(emb_outputs[2])
    dist_out = cluster_distance_calculations(calc_dist_cluster_input)
    emit:
        distances = dist_out

}
workflow{
    println "current directory ${projectDir}"
    println "Output dir:  ${params.output_base}"
    println "DNABERT-S dir:     ${params.test_model_dir}"
    data_dir = "${projectDir.parent}/data"

    input_ch = channel.fromPath(params.input_txt)

    fasta_files = download_fasta(input_ch)
    processed_fastas = process_fasta(fasta_files.collect(), input_ch)
    // Calculate coverage
    coverage_file_SHORT = calculate_coverage_for_short_reads(processed_fastas[0].collect(), input_ch)

    short_ch = short_reads(processed_fastas[1], input_ch, coverage_file_SHORT)
    long_ch  = long_reads(processed_fastas[1], input_ch)

    reads_ch = short_ch.mix(long_ch)
    distances_ch = process_read_pipeline(reads_ch).distances
    // distances_ch: tuple(read_type, npz)

    short_npz = distances_ch.filter{ it[0] == 'short_reads' }.map{ it[1] }.first()
    long_npz  = distances_ch.filter{ it[0] == 'long_reads'  }.map{ it[1] }.first()
    cluster_distance_visualization(short_npz,long_npz)
}