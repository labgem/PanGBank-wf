process PARSE_GENOMES_AND_TAXONOMY {
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    path genomes
    path taxonomy
    val min_genomes

    output:
    // path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    // when:
    // task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in labgem/pangbank/bin/
    """
    parse_genomes_and_taxonomy.py --genomes $genomes --taxonomy $taxonomy\
                                    --min_genomes $min_genomes --species_summary_file species_summary.tsv\
                                    --ppanggolin_files_outdir ppanggolin_input_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

