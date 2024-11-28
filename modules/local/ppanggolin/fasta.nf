process PPANGGOLIN_FASTA {
    tag "${meta.species}"

    cpus 1

    memory "2GB"

    conda "bioconda::ppanggolin>=2.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.0.4--py310h4b81fae_0' :
        'biocontainers/ppanggolin:2.0.4--py310h4b81fae_0' }"

    input:
        tuple val(meta), path(pangenome)

    output:
        tuple val(meta), path("${meta.species}/persistent_nucleotide_families.fasta.gz") , emit: persistent_families_fasta
        tuple val(meta), path("${meta.species}/all_protein_families.faa.gz") , emit: all_families_faa
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    ppanggolin fasta -p $pangenome  -o ${meta.species} --gene_families persistent --prot_families all --compress

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ppanggolin: \$(ppanggolin --version | sed 's/ppanggolin //g')
    END_VERSIONS
    """

}

