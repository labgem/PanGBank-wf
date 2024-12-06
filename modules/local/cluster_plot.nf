process CLUSTER_PLOT {
    label 'process_low'

    conda "bioconda::ppanggolin>=2.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.0.4--py310h4b81fae_0' :
        'biocontainers/ppanggolin:2.0.4--py310h4b81fae_0' }"

    input:
    path cluster_stat
    path distance_count

    output:
    path "distance_plots/*.html"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_cluster_stat.py \\
            --cluster_stat $cluster_stat \\
            --distance_to_count $distance_count \\
            --output_dir distance_plots
    """
}
