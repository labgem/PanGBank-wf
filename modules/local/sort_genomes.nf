process SORT_GENOMES {
    label 'process_single'
    tag "${meta.id}"

    // reuse ppanggolin env as it as already been downloaded and used
    conda "bioconda::ppanggolin=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.2.3--hbcbf7aa_0' :
        'biocontainers/ppanggolin:2.2.3--hbcbf7aa_0' }"

    input:
        tuple val(meta), path(genome_stat_file)

    output:
        tuple  val(meta), path("sorted_genomes.txt.gz")       , emit: sorted_genomes_list
        path "sorted_genomes.txt.gz"                          , emit: sorted_genomes_stats
        path "versions.yml"      , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:

    """
    sort_genomes.py --genome_stats $genome_stat_file \\
                    --sort_by L90 L75 L50 auN \\
                    --sorted_genome_list sorted_genomes.txt.gz \\
                    --sorted_genome_stats sorted_genomes.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """

}
