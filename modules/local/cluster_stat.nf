process CLUSTER_STAT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::ppanggolin>=2.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.0.4--py310h4b81fae_0' :
        'biocontainers/ppanggolin:2.0.4--py310h4b81fae_0' }"

    input:
    tuple val(meta), path(cluster_composition), path(phylip_matrix)

    output:
    path "cluster_stat.tsv", emit: cluster_stat
    path "distance_count.tsv", emit: distance_count

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    compute_cluster_stat.py \\
            --cluster_composition $cluster_composition \\
            --phylip_matrix $phylip_matrix \\
            --species $meta.id \\
            --cluster_stat cluster_stat.tsv \\
            --distance_count_file distance_count.tsv
    """
}
