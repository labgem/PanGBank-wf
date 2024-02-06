process PPANGGOLIN {
    label 'process_single'
    tag "${meta.species} - (${meta.genomes_count})"

    conda "bioconda::ppanggolin>=2.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.0.2--py39hf95cd2a_0' :
        'biocontainers/ppanggolin:2.0.2--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(genome_file), path("genomes/*")

    output:
    // path 'ppanggolin_input_files/*.tsv'       , emit: ppanggo_inputs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ppanggolin all --anno $genome_file --output ppanggolin_results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ppanggolin: \$(ppanggolin --version | sed 's/ppanggolin //g')
    END_VERSIONS
    """
}

