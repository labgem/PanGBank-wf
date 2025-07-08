process PARSE_GENOMES_AND_TAXONOMY {
    label 'process_single'

    // reuse ppanggolin env as it as already been downloaded and used
    conda "bioconda::ppanggolin=2.2.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin:2.2.4--h0fa9677_0' :
        'biocontainers/ppanggolin:2.2.4--h0fa9677_0' }"

    input:
    path genomes
    path taxonomy
    path genome_metadata
    val min_genomes

    output:
    path "ppanggolin_input_files/*/input_genomes.tsv.gz"       , emit: ppanggo_inputs
    path "species_summary.tsv"                                 , emit: summary
    path "ppanggolin_input_files/*/genomes_metadata.tsv.gz"    , optional: true, emit: genome_metadata
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome_metadata_arg = genome_metadata.name != 'NO_FILE' ? "--genome_metadata $genome_metadata" : ''

    """
    parse_genomes_and_taxonomy.py --genomes $genomes --taxonomy $taxonomy\
                                    --min_genomes $min_genomes --species_summary_file species_summary.tsv\
                                    --outdir ppanggolin_input_files $genome_metadata_arg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
