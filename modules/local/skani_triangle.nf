process SKANI_TRIANGLE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::skani=0.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/skani:0.3.1--ha6fb395_0':
        'biocontainers/skani:0.3.1--ha6fb395_0' }"

    input:
    tuple val(meta), path(queries), path(sorted_genomes)

    output:
    tuple val(meta), path("distance.phylip") , emit: phylip_matrix
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    skani \\
        triangle \\
            -l ${queries} \\
            --distance \\
            -t ${task.cpus} \\
            ${args} \\
        | awk 'FNR==NR {idx[\$0]=FNR-1; next} FNR==1 {print; next} {\$1=idx[\$1]; print}' \\
            <(zcat ${sorted_genomes}) - \\
        > distance.phylip
    # The awk command converts genome names to indices based on the order in the sorted_genomes file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version 2>&1 | sed 's/^.*skani //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch distance.phylip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version 2>&1 | sed 's/^.*skani //; s/ .*\$//')
    END_VERSIONS
    """
}
