process MD5SUM_ON_FILES {

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
        tuple val(meta), path(genome_file)

    output:
        tuple  val(meta), path("genomes_md5sum.tsv")       , emit: genomes_md5sum_files
        path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Assign input and output file paths
    output_file=genomes_md5sum.tsv

    # Create or clear the output file
    echo -e "name\tfile_name\tmd5_sum" > "\$output_file"

    # Read the input file line by line
    while IFS=\$'\t' read -r genome_name genome_path; do


        # Calculate the MD5 checksum of the genome file
        md5_sum=\$(md5sum "\$genome_path" | awk '{ print \$1 }')

        # Get the name of the genome file
        genome_filename=\$(basename "\$genome_path")

        # Write the output to the output file
        echo -e "\${genome_name}\t\${genome_filename}\t\${md5_sum}" >> "\$output_file"

    done < "$genome_file"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        md5sum: \$(echo \$(md5sum --version 2>&1 | head -n 1| sed 's/^.*) //;' ))
    END_VERSIONS

    """
}
