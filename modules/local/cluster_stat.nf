process CLUSTER_STAT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::ppanggolin=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin:2.2.1--py311haab0aaa_1' :
        'biocontainers/ppanggolin:2.2.1--py311haab0aaa_1' }"

    input:
    tuple val(meta), path(cluster_composition), path(phylip_matrix)

    output:
    tuple val(meta), path("cluster_stat.tsv"), emit: cluster_stat
    tuple val(meta), path("distance_count.tsv"), emit: distance_count
    path "versions.yml"      , emit: versions


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



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS

    """
}
