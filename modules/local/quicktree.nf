process QUICKTREE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::quicktree"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quicktree%3A2.2--h470a237_1':
        'biocontainers/quicktree%3A2.2--h470a237_1' }"


    input:
    tuple val(meta), path(phylip_matrix)

    output:
    tuple val(meta), path("tree.nw")    ,  emit: tree
    path "versions.yml"               ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    quicktree -in m $phylip_matrix > tree.nw


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quicktree: \$(quicktree -v | sed 's/quicktree //g')
    END_VERSIONS
    """

}
