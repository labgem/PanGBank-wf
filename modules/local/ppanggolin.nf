process PPANGGOLIN {
    // label 'process_single'
    tag "${meta.species} - (${meta.genomes_count})"

    // A dynamic label would be perfect here but does not work.. https://github.com/nextflow-io/nextflow/issues/894

    queue { meta.genomes_count > params.large_pangenome_cutoff ? params.large_pangenome_queue : params.regular_pangenome_queue }

    time { meta.genomes_count > params.large_pangenome_cutoff ? '5days' : '23:50:00' }
    // clusterOptions { meta.genomes_count > 5000 ? '--tmp 50G --exclusive=user' : ''  } // node with at least XGo and exclusif to the user

    // 16 cpu when more than 5k, from 1 to 16cpu from 1 to 5k genomes
    cpus { meta.genomes_count > params.large_pangenome_cutoff ? "16" : "${Math.round(Math.ceil(meta.genomes_count / 312))}" }

    // With >5K  genomes : 30GB per cpu otherwise 8GB/cpu
    memory { meta.genomes_count > params.large_pangenome_cutoff ?  "${16*30}GB" : "${Math.ceil(2 + (meta.genomes_count / 312)*8)}GB" }
    // memory { meta.genomes_count > 30 ? '1 GB' : '3 GB' }

    tag { "${meta.species} ${meta.genomes_count}genomes ${Math.ceil((meta.genomes_count / 312)*8)}GB ${Math.round(Math.ceil(meta.genomes_count / 312))}cpus" }

    conda "bioconda::ppanggolin>=2.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.0.3--py39hf95cd2a_0' :
        'biocontainers/ppanggolin:2.0.3--py310h4b81fae_0' }"

    input:
    tuple val(meta), path(genome_file), path("genomes/*")


    output:
    tuple  val(meta), path("${meta.species}/pangenome.h5")       , emit: pangenome
    // tuple  val(meta), path("${meta.species}/info.yaml")          , emit: pangenome_info

    path "${meta.species}.yaml"                                  , emit: pangenome_info
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
    ppanggolin all $input --output ${meta.species} --no_flat_files  --cpu $task.cpus  $tmpdir

    ppanggolin info -p ${meta.species}/pangenome.h5 --content > ${meta.species}.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ppanggolin: \$(ppanggolin --version | sed 's/ppanggolin //g')
    END_VERSIONS
    """

}

