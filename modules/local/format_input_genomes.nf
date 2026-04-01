process FORMAT_INPUT_GENOMES {
    label 'process_single'
    tag "$meta.id"

    // reuse ppanggolin env as it as already been downloaded and used
    conda "bioconda::ppanggolin=2.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.3.0--py312h247cb63_0' :
        'biocontainers/ppanggolin:2.3.0--py312h247cb63_0' }"

    input:
    tuple val(meta), path(selected_genomes), path(genome_name_to_fasta), path(fasta_to_original_path)
    path reference_genomes

    output:
    tuple val(meta), path("formatted_selected_genomes.tsv.gz"), emit: genome_selection_ppanggo_input
    path "summary_selection.tsv", emit: dereplication_summary
    path "versions.yml"      , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def fasta_to_original_path_arg = fasta_to_original_path.name != 'NO_FILE' ? "--fasta_to_original_path $fasta_to_original_path" : ''
    def reference_genomes_arg = reference_genomes ? "--reference_genomes $reference_genomes" : ''

    """
    echo $reference_genomes
    format_ppanggo_input_from_genome_selection.py --selected_genomes $selected_genomes \\
                                --genome_name_to_path $genome_name_to_fasta \\
                                --formatted_genomes formatted_selected_genomes.tsv.gz \\
                                $fasta_to_original_path_arg $reference_genomes_arg \\
                                --summary_selection summary_selection.tsv --species ${meta.id} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

}
