process MASH_DIST_TO_PHYLIP {
    tag "$meta.id"

    // search a mulled container with mash and python as describe here: https://nf-co.re/docs/guidelines/components/modules#re-use-of-multi-tool-containers
    // mulled-search  -d quay singularity  --search mash python | grep mulled
    conda "bioconda::mash=2.3 anaconda::numpy=2.1.3 conda-forge::tqdm"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-83c3b71f1ea527155742db36d554b890aaa6a3d4:06958a0680d045f715a4f8b3bcce64df31d13ef8-0' :
        'biocontainers/mulled-v2-83c3b71f1ea527155742db36d554b890aaa6a3d4:06958a0680d045f715a4f8b3bcce64df31d13ef8-0' }"

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
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
