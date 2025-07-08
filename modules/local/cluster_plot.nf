process CLUSTER_PLOT {
    label 'process_single'

    conda "bioconda::ppanggolin=2.2.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin:2.2.4--h0fa9677_0' :
        'biocontainers/ppanggolin:2.2.4--h0fa9677_0' }"

    input:
    path cluster_stat
    path distance_count

    output:
    path "plots/*.html", emit: plots
    path "versions.yml"      , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_cluster_stat.py \\
            --cluster_stat $cluster_stat \\
            --distance_to_count $distance_count \\
            --output_dir plots

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        plotly: \$(python -c "import plotly; print(plotly.__version__)")
        scipy: \$(python -c "import scipy; print(scipy.__version__)")
    END_VERSIONS
    """
}
