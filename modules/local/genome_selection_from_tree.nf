process GENOME_SELECTION_FROM_TREE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::treeswift=1.1.45"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/treeswift:1.1.45--pyh7e72e81_0':
        'biocontainers/treeswift:1.1.45--pyh7e72e81_0' }"


    input:
    tuple val(meta), path(tree), path(sorted_genomes)
    val number_of_genomes

    output:
    tuple val(meta), path("selected_genomes.txt")    ,  emit: selected_genomes
    tuple val(meta), path("cluster_composition.txt")    ,  emit: cluster_composition
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    select_genomes_from_tree.py --tree $tree \\
                                --number_of_genomes $number_of_genomes \\
                                --sorted_genomes $sorted_genomes \\
                                --cluster_composition "cluster_composition.txt" \\
                                --selected_genomes selected_genomes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        treeswift: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('treeswift').version)")
    END_VERSIONS
    """

}
