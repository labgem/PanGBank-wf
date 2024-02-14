process PPANGGOLIN {
    // label 'process_single'
    tag "${meta.species} - (${meta.genomes_count})"

    // A dynamic label would be perfect here but does not work.. https://github.com/nextflow-io/nextflow/issues/894

    conda "bioconda::ppanggolin>=2.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.0.2--py39hf95cd2a_0' :
        'biocontainers/ppanggolin:2.0.2--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(genome_file), path("genomes/*")


    output:
    // path 'ppanggolin_input_files/*.tsv'       , emit: ppanggo_inputs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def input = meta.file_type == "annotation" ? "--anno $genome_file" : "--fasta $genome_file"
    """
    ppanggolin all $input --output ppanggolin_results --no_flat_files  --cpu $task.cpus --tmpdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ppanggolin: \$(ppanggolin --version | sed 's/ppanggolin //g')
    END_VERSIONS
    """
}

