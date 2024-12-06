process SORT_GENOMES {
    // label 'process_single'
    tag "${meta.id}"

    // reuse ppanggolin env as it as already been downloaded and used
    conda "bioconda::ppanggolin=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.2.1--py311haab0aaa_1' :
        'biocontainers/ppanggolin%3A2.2.1--py311haab0aaa_1' }"

    input:
        tuple val(meta), path(genome_stat_file)

    output:
        tuple  val(meta), path("sorted_genomes.txt")       , emit: sorted_genomes_list
        path "sorted_genomes.tsv"                          , emit: sorted_genomes_stats

    when:
        task.ext.when == null || task.ext.when

    script:

    """
    sort_genomes.py --genome_stats $genome_stat_file \\
                    --sort_by L90 L75 L50 auN \\
                    --sorted_genome_list sorted_genomes.txt \\
                    --sorted_genome_stats sorted_genomes.tsv

    """

}

