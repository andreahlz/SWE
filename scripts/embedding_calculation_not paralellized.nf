#!/usr/bin/env nextflow

/*
* Apply Insilicoseq to fasta file
*/ 
process iss {
    publishDir "${projectDir.parent}/data/fastq_files", mode: 'copy'

    input:
        file fasta_file
    output:
        path "${fasta_file.baseName}*"
    script:
        """
        echo 'starting read generation for ${fasta_file.baseName}'
        echo 'fasta_file name = ${fasta_file}'
        iss generate \
            --genomes ${fasta_file} \
            --model MiSeq \
            --n_reads 2000000 \
            --output ${fasta_file.baseName}
        """
}
process fastq_to_tsv{
    publishDir "${projectDir.parent}/data/tsvs",mode:'copy'

    input:
        path fastq_file
        path py_script
        val fastq_dir
    output:
        path "${fastq_file.baseName.replace("_genomic_R1","")}.tsv"
        //path "${fastq_file.baseName}.tsv"
    script:
        """
        python3 ${py_script} "${fastq_dir}/${fastq_file}"
        """
}
process combine_tsvs{
    publishDir "${projectDir.parent}/data/DNABERTS_input/gut_microbiome_0",mode:'copy'

    input:
        path tsv_files
    output:
        file "DNABERTS_input/gut_microbiome_0/all_tsvs.tsv"
    script: 
    """
    echo "$tsv_files"
    mkdir DNABERTS_input/gut_microbiome_0 -p
    cat ${tsv_files} > "DNABERTS_input/gut_microbiome_0/all_tsvs.tsv"
    """
}
process DNABERTS{
    publishDir "${projectDir.parent}/data/DNABERTS_results",mode:'copy'

    input:
        path py_file
        path test_model_dir
        path data_dir
    output:
        file "embeddings_visualization.png"
        path "data/embeddings/gut_microbiome_0/clustering_all_tsvs/test.npy"
    script:
    """
    export TRANSFORMERS_VERBOSITY=error
    python3 ${py_file} \
        --test_model_dir=${test_model_dir} \
        --model_list "test" \
        --data_dir ${data_dir}
    """
}
workflow {
    println "Script directory: ${baseDir}"
    def parentDir = projectDir.parent
    // channel containing all fasta files
    fasta_files = channel.fromPath("${parentDir}/data/genomes/processed/*")
    iss(fasta_files)
    // Python-Skript als Input
    py_script = file("${projectDir}/fastq_to_tsv_local.py")
    
    // Channel mit allen R1 FASTQ-Files
    //fastq_files = channel.fromPath("${projectDir}/data/fastq_files/*_R1.fastq")
    
    // Prozess ausführen
    //fastq_to_tsv(fastq_files, py_script, "${parentDir}/data/fastq_files")
    //def alltsvs = fastq_to_tsv.out.collect()
    processed_fastq = iss.out.flatten().filter(~ /.*_R1\.fastq/)
    fastq_to_tsv(processed_fastq, py_script, "${parentDir}/data/fastq_files")
    combine_tsvs(fastq_to_tsv.out.collect())
    py_file = file("${projectDir}/eval_clustering_classification_changed.py")
    println "${parentDir}/data/DNABERTS_input"
    test = combine_tsvs.out.map{ file -> file.parent.parent }
    test.view()
    DNABERTS(py_file,"/home/barbara/evaluate/DNABERT-S",test)
    //DNABERTS(py_file,"/home/barbara/evaluate/DNABERT-S","${parentDir}/data/DNABERTS_input")
}