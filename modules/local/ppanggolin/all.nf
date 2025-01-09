process PPANGGOLIN_ALL {
    // label 'process_single'
    tag "${meta.species} - (${meta.genomes_count})"

    // A dynamic label would be perfect here but does not work.. https://github.com/nextflow-io/nextflow/issues/894

    queue { meta.genomes_count > params.large_pangenome_cutoff ? params.large_pangenome_queue : params.regular_pangenome_queue }

    time { meta.genomes_count >= params.large_pangenome_cutoff ?
            '23:50:00' :
            '12:00:00'}
    // clusterOptions { meta.genomes_count > 5000 ? '--tmp 50G --exclusive=user' : ''  } // node with at least XGo and exclusif to the user

    // 16 cpu when more than 5k, from 1 to 16cpu from 1 to 5k genomes
    cpus {
        meta.genomes_count >= params.large_pangenome_cutoff ?
        "16" :
        (meta.genomes_count >= (params.large_pangenome_cutoff / 2)) ?
        "8" :
        "4"
    }

    memory {
        meta.genomes_count >= params.large_pangenome_cutoff ?
        "${8*16}GB" :
        (meta.genomes_count >= (params.large_pangenome_cutoff / 2)) ?
        "${8*8}GB" :
        "${4*4}GB"
    }

    tag { "${meta.species} ${meta.genomes_count}genomes ${Math.ceil((meta.genomes_count / 312)*8)}GB ${Math.round(Math.ceil(meta.genomes_count / 312))}cpus" }

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
        path "${meta.species}/genomes_statistics.tsv*"

        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def input = meta.file_type == "annotation" ? "--anno $genome_file" : "--fasta $genome_file"

    def tmpdir = ""
    if (params.large_pangenome_tmpdir && meta.genomes_count > params.large_pangenome_cutoff){
        tmpdir =  "--tmpdir ${params.large_pangenome_tmpdir}"
    }
    else if (params.regular_pangenome_tmpdir && meta.genomes_count <= params.large_pangenome_cutoff) {
        tmpdir =  "--tmpdir ${params.regular_pangenome_tmpdir}"
    }

    """
    ppanggolin all $input --output ${meta.species} --config $ppanggolin_config  --cpu $task.cpus  $tmpdir

    ppanggolin info -p ${meta.species}/pangenome.h5 --content > ${meta.species}.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ppanggolin: \$(ppanggolin --version | sed 's/ppanggolin //g')
    END_VERSIONS
    """

}

