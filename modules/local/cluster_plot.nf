process CLUSTER_PLOT {
    label 'process_low'

    conda "bioconda::ppanggolin=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin:2.2.1--py311haab0aaa_1' :
        'biocontainers/ppanggolin:2.2.1--py311haab0aaa_1' }"

    input:
    path cluster_stat
    path distance_count

    output:
    path "plots/*.html"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_cluster_stat.py \\
            --cluster_stat $cluster_stat \\
            --distance_to_count $distance_count \\
            --output_dir plots
    """
}
