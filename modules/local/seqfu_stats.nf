process SEQFU_STATS_FROM_FILE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::seqfu=1.22.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqfu:1.22.3--h1eb128b_0':
        'biocontainers/seqfu:1.22.3--h1eb128b_0' }"


    input:
    // stats can get one or more fasta or fastq files
    tuple val(meta), path(genome_path_file)

    output:
    tuple val(meta), path("genome_stat.tsv")    ,  emit: stats
    path "versions.yml"               ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    # Clear the output file or create it fresh
    > "genome_stat.tsv"

    # Loop over each genome path in the list
    while IFS= read -r genome_path; do

        arg="--noheader"
        # If the output file is empty, include the header, otherwise append with --noheader
        if [ ! -s "genome_stat.tsv" ]; then
            arg=" "
        fi

        seqfu stats --index \$arg "\$genome_path" >> "genome_stat.tsv"

    done < $genome_path_file


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqfu: \$(seqfu version)
    END_VERSIONS
    """

}
