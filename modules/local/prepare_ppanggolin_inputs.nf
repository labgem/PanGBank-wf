process PREPARE_PPANGGOLIN_INPUTS {
    label 'process_single'

    // reuse ppanggolin env as it as already been downloaded and used
    conda "bioconda::ppanggolin=2.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.3.0--py312h247cb63_0' :
        'biocontainers/ppanggolin:2.3.0--py312h247cb63_0' }"

    input:
    path genomes
    path taxonomy
    path genome_metadata
    path translation_tables
    path gtdb_cluster
    val min_genomes

    output:
    path "ppanggolin_input_files/*/input_genomes.tsv.gz"       , emit: ppanggo_inputs
    path "species_summary.tsv"                                 , emit: summary
    path "ppanggolin_input_files/*/genomes_metadata.tsv.gz"    , optional: true, emit: genome_metadata
    path "ppanggolin_input_files/species_to_translation_tables.tsv"    , optional: true, emit: species_translation_tables
    path "ppanggolin_input_files/*/pangenome_taxonomy.txt"
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome_metadata_arg = genome_metadata.name != 'NO_FILE' ? "--genome_metadata $genome_metadata" : ''
    def genome_translation_tables_arg = translation_tables.name != 'NO_FILE_2' ? "--genome_translation_table $translation_tables" : ''
    def gtdb_cluster_arg = gtdb_cluster.name != 'NO_FILE_3' ? "--species_to_merge $gtdb_cluster" : ''

    """
    prepare_ppanggolin_inputs.py --genomes $genomes --taxonomy $taxonomy\
                                    --min_genomes $min_genomes --species_summary_file species_summary.tsv\
                                    --outdir ppanggolin_input_files \
                                    $genome_metadata_arg $genome_translation_tables_arg $gtdb_cluster_arg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
