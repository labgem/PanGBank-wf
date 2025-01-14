process PARSE_GENOMES_AND_TAXONOMY {
    label 'process_single'

    // reuse ppanggolin env as it as already been downloaded and used
    conda "bioconda::ppanggolin=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin:2.2.1--py311haab0aaa_1' :
        'biocontainers/ppanggolin:2.2.1--py311haab0aaa_1' }"
    input:
    path genomes
    path taxonomy
    val min_genomes

    output:
    path "ppanggolin_input_files/*.tsv"       , emit: ppanggo_inputs
    path "species_summary.tsv"                , emit: summary
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in labgem/pangbank/bin/
    """
    parse_genomes_and_taxonomy.py --genomes $genomes --taxonomy $taxonomy\
                                    --min_genomes $min_genomes --species_summary_file species_summary.tsv\
                                    --ppanggolin_files_outdir ppanggolin_input_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}

