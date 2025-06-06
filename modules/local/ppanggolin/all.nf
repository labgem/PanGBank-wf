process PPANGGOLIN_ALL {
    label 'process_medium'


    tag { "${meta.species}" }

    conda "bioconda::ppanggolin=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.2.3--hbcbf7aa_0' :
        'biocontainers/ppanggolin:2.2.3--hbcbf7aa_0' }"

    input:
        tuple val(meta), path(genome_file)
        path(ppanggolin_config)

    output:
        tuple  val(meta), path("${meta.species}/pangenome.h5")          , emit: pangenome
        tuple  val(meta), path("${meta.species}/genomes_statistics.tsv.gz")  , emit: genomes_statistics
        path "${meta.species}.yaml"                                     , emit: pangenome_info
        path "${meta.species}/metadata/*.tsv.gz"                        , optional: true
        path "${meta.species}/proksee/*.json*"
        path "versions.yml"                                             , emit: versions
        path "${meta.species}/*.html"

    when:
        task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_arg = meta.file_type == "annotation" ? "--anno" : "--fasta"


    """
    ppanggolin all $input_arg  $genome_file --output ${meta.species} --config $ppanggolin_config  --cpu $task.cpus  $args

    ppanggolin info --pangenome ${meta.species}/pangenome.h5 --content > ${meta.species}.yaml

    # write metadata collected on annotation if any
    ppanggolin write_metadata --pangenome ${meta.species}/pangenome.h5 --output ${meta.species}/metadata/ --compress

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ppanggolin: \$(ppanggolin --version | sed 's/ppanggolin //g')
    END_VERSIONS
    """

}

