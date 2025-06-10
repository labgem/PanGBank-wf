process PPANGGOLIN_FASTA {
    tag "${meta.species}"
    label 'process_single'

    conda "bioconda::ppanggolin=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.2.3--hbcbf7aa_0' :
        'biocontainers/ppanggolin:2.2.3--hbcbf7aa_0' }"

    input:
        tuple val(meta), path(pangenome)

    output:
        // tuple val(meta), path("${meta.species}/persistent_nucleotide_families.fasta.gz") , emit: persistent_families_fasta
        tuple val(meta), path("persistent_nucleotide_families/${meta.species}.fasta.gz") , emit: persistent_families_fasta
        tuple val(meta), path("${meta.species}/all_protein_families.faa.gz") , emit: all_families_faa
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    ppanggolin fasta -p $pangenome  -o ${meta.species} --gene_families persistent --prot_families all --compress
    mkdir -p persistent_nucleotide_families
    mv ${meta.species}/persistent_nucleotide_families.fasta.gz persistent_nucleotide_families/${meta.species}.fasta.gz


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ppanggolin: \$(ppanggolin --version | sed 's/ppanggolin //g')
    END_VERSIONS
    """

}
