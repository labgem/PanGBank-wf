process PARSE_GENOMES_AND_TAXONOMY {
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/python:3.10' :
    //     'biocontainers/python:3.10' }"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.0.4--py310h4b81fae_0' :
        'biocontainers/ppanggolin:2.0.4--py310h4b81fae_0' }"
    input:
    path genomes
    path taxonomy
    val min_genomes

    output:
    path "ppanggolin_input_files/*.tsv"       , emit: ppanggo_inputs
    path "species_summary.tsv"

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in labgem/pangbank/bin/
    """
    parse_genomes_and_taxonomy.py --genomes $genomes --taxonomy $taxonomy\
                                    --min_genomes $min_genomes --species_summary_file species_summary.tsv\
                                    --ppanggolin_files_outdir ppanggolin_input_files

    """
}

