// Publish the final input genomes file to the results directory.
// This is a dedicated publish step so that there is a single, deterministic
// source of truth for pangenomes/<species>/input_genomes.tsv.gz —
// regardless of whether the species went through dereplication or not.
process PUBLISH_INPUT_GENOMES {
    tag "$meta.species"
    label 'process_single'

    conda "bioconda::ppanggolin=2.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.3.0--py312h247cb63_0' :
        'biocontainers/ppanggolin:2.3.0--py312h247cb63_0' }"

    input:
    tuple val(meta), path(input_genomes)

    output:
    path "input_genomes.tsv.gz", emit: input_genomes

    script:
    """
    cp -L $input_genomes input_genomes.tsv.gz

    """
}
