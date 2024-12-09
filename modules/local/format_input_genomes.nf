process FORMAT_INPUT_GENOMES {
    label 'process_single'
    tag "$meta.id"

    // reuse ppanggolin env as it as already been downloaded and used
    conda "bioconda::ppanggolin=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin:2.2.1--py311haab0aaa_1' :
        'biocontainers/ppanggolin:2.2.1--py311haab0aaa_1' }"

    input:
    tuple val(meta), path(selected_genomes), path(genome_name_to_fasta), path(fasta_to_original_path)
    path reference_genomes

    output:
    tuple val(meta), path("formatted_selected_genomes.tsv"), emit: genome_selection_ppanggo_input

    when:
    task.ext.when == null || task.ext.when

    script:
    def fasta_to_original_path_arg = fasta_to_original_path.name != 'NO_FILE' ? "--fasta_to_original_path $fasta_to_original_path" : ''
    def reference_genomes_arg = reference_genomes.name ? "--reference_genomes $reference_genomes" : ''

    """
    echo $reference_genomes
    format_ppanggo_input_from_genome_selection.py --selected_genomes $selected_genomes \\
                                --genome_name_to_path $genome_name_to_fasta \\
                                --formatted_genomes formatted_selected_genomes.tsv \\
                                $fasta_to_original_path_arg $reference_genomes_arg

    """

}

