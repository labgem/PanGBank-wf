process PPANGGOLIN_ALL {
    label 'process_medium'


    tag { "${meta.species}" }

    conda "bioconda::ppanggolin=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin:2.2.1--py311haab0aaa_1' :
        'biocontainers/ppanggolin:2.2.1--py311haab0aaa_1' }"

    input:
        tuple val(meta), path(genome_file)
        path(ppanggolin_config)

    output:
        tuple  val(meta), path("${meta.species}/pangenome.h5")       , emit: pangenome
        path "${meta.species}.yaml"                                  , emit: pangenome_info
        path "${meta.species}/genomes_statistics.tsv.gz"             , emit: genomes_statistics
        path "${meta.species}/*.html"

        path "${meta.species}/metadata/*.tsv.gz"                     , optional: true

        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def input_arg = meta.file_type == "annotation" ? "--anno" : "--fasta"

    def tmpdir = ""
    if (params.large_pangenome_tmpdir && meta.genomes_count > params.large_pangenome_cutoff){
        tmpdir =  "--tmpdir ${params.large_pangenome_tmpdir}"
    }
    else if (params.regular_pangenome_tmpdir && meta.genomes_count <= params.large_pangenome_cutoff) {
        tmpdir =  "--tmpdir ${params.regular_pangenome_tmpdir}"
    }

    // TODO remove zcat when genome input will be supported by ppanggolin

    """
    zcat $genome_file > genomes_file_list.txt
    ppanggolin all $input_arg  genomes_file_list.txt --output ${meta.species} --config $ppanggolin_config  --cpu $task.cpus  $tmpdir

    ppanggolin info --pangenome ${meta.species}/pangenome.h5 --content > ${meta.species}.yaml

    # write metadata collected on annotation if any
    ppanggolin write_metadata --pangenome ${meta.species}/pangenome.h5 --output ${meta.species}/metadata/ --compress

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ppanggolin: \$(ppanggolin --version | sed 's/ppanggolin //g')
    END_VERSIONS
    """

}

