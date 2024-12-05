process ANY2FASTA {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::any2fasta=0.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/any2fasta:0.4.2--hdfd78af_3':
        'biocontainers/any2fasta:0.4.2--hdfd78af_3' }"


    input:
    tuple val(meta), path(genome_path_file)

    output:
    tuple val(meta), path("genome_name_to_fasta.tsv")    ,  emit: genome_path_fasta
    tuple val(meta), path("fasta_to_orginal_path.tsv")    ,  emit: fasta_to_orginal_path
    path "versions.yml"               ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    while IFS=\$'\t' read -r genome_name genome_path; do
        # Run the any2fasta command and compress the output
        any2fasta -q -u "\$genome_path" | gzip > "\${genome_name}.fasta.gz"
        fasta_genome_path=`realpath "\${genome_name}.fasta.gz"`

        echo -e \$genome_name'\\t'\$fasta_genome_path >> genome_name_to_fasta.tsv
        echo -e \$fasta_genome_path'\\t'\$genome_path >> fasta_to_orginal_path.tsv


    done < $genome_path_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        any2fasta: \$(any2fasta -v | sed 's/any2fasta //g')
    END_VERSIONS
    """

}
