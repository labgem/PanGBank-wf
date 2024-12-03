process MASH_DIST_TO_PHYLIP {
    tag "$meta.id"
    label 'process_low'

    // search a mulled container with mash and python as describe here: https://nf-co.re/docs/guidelines/components/modules#re-use-of-multi-tool-containers
    // mulled-search  -d quay singularity  --search mash python | grep mulled
    conda "bioconda::mash=2.3 conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        ' https://depot.galaxyproject.org/singularity/mulled-v2-7f933c366aa077c57320ffc29bf8f7e2d0ea9234:f2c964aa6075b63b60450d37ba79b2f1fb0f6020-0' :
        'biocontainers/mulled-v2-7f933c366aa077c57320ffc29bf8f7e2d0ea9234:f2c964aa6075b63b60450d37ba79b2f1fb0f6020-0' }"

    input:
    tuple val(meta), path(mash_sketch), path(sorted_genomes)

    output:
    tuple val(meta), path("distance.phylip"), emit: phylip_matrix
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mash \\
        dist \\
        -p $task.cpus \\
        $mash_sketch \\
        $mash_sketch | \\
        mash_dist_to_phylip_matrix.py \\
                --sorted_genomes_file $sorted_genomes \\
                --mash_dist_result - \\
                --phylip_matrix distance.phylip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """
}
