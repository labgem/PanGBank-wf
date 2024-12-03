process GENOME_SELECTION_FROM_TREE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::treeswift"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/treeswift%3A1.1.0--py_0':
        'biocontainers/treeswift%3A1.1.0--py_0' }"


    input:
    tuple val(meta), path(tree), path(sorted_genomes),  path(genome_name_to_path)
    val number_of_genomes

    output:
    tuple val(meta), path("selected_genomes.tsv")    ,  emit: selected_genomes
    tuple val(meta), path("cluster_composition.txt")    ,  emit: cluster_composition


    when:
    task.ext.when == null || task.ext.when

    script:
    """

    select_genomes_from_tree.py --tree $tree \\
                                --number_of_genomes $number_of_genomes \\
                                --sorted_genomes $sorted_genomes \\
                                --genome_name_to_path $genome_name_to_path \\
                                --cluster_composition "cluster_composition.txt" \\
                                --selected_genomes selected_genomes.tsv

    """

}
