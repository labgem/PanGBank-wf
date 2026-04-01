process PPANGGOLIN_METADATA {
    label 'process_medium'


    tag { "${meta.species}" }

    conda "bioconda::ppanggolin=2.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.3.0--py312h247cb63_0' :
        'biocontainers/ppanggolin:2.3.0--py312h247cb63_0' }"

    input:
    tuple val(meta), path(pangenome), path(genome_metadata)

    output:
    path "versions.yml", emit: versions
    path "pangenome_with_metdata.h5", emit: pangenome_with_metadata

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    cp ${pangenome} pangenome_with_metdata.h5 # required to copy to not change the original file to avoid issues when resuming

    ppanggolin metadata -p pangenome_with_metdata.h5 --metadata ${genome_metadata} \\
                        --source "pangbank_wf"  --assign genomes \\
                        --omit --force  # force in case of resume and omit in case genome are missing because of dereplication



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ppanggolin: \$(ppanggolin --version | sed 's/ppanggolin //g')
    END_VERSIONS
    """
}
