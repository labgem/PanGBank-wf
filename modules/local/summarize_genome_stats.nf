process SUMMARIZE_GENOME_STATS {
    label 'process_single'
    tag "${meta.species}"


    // reuse ppanggolin env as it as already been downloaded and used
    conda "bioconda::ppanggolin=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.2.3--hbcbf7aa_0' :
        'biocontainers/ppanggolin:2.2.3--hbcbf7aa_0' }"

    input:
    tuple val(meta), path(genome_stat_file)

    output:
    path  "${meta.species}.yaml", emit: genome_stats_summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    summarize_genome_stats.py --genome_stats $genome_stat_file --output ${meta.species}.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}

