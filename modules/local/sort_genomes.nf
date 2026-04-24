process SORT_GENOMES {
    label 'process_single'
    tag "${meta.id}"

    // reuse ppanggolin env as it as already been downloaded and used
    conda "bioconda::ppanggolin=2.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.3.0--py312h247cb63_0' :
        'biocontainers/ppanggolin:2.3.0--py312h247cb63_0' }"

    input:
        tuple val(meta), path(genome_name_to_path), path(genome_stat_file), path(genome_metadata_file)

    output:
        tuple  val(meta), path("sorted_genomes.txt.gz")       , emit: sorted_genomes_list
        path "sorted_genomes.tsv.gz"                          , emit: sorted_genomes_stats
        path "versions.yml"      , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:

    def genome_metadata_args = genome_metadata_file.name != 'NO_FILE' ? "--genome_metadata $genome_metadata_file --completeness_sorting" : ''

    """
    sort_genomes.py --genome_name_to_path $genome_name_to_path \\
                    --genome_stats $genome_stat_file \\
                    --sorted_genome_list sorted_genomes.txt.gz \\
                    --sorted_genome_stats sorted_genomes.tsv.gz \\
                    $genome_metadata_args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

}
