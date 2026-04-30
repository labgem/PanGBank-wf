process FIND_GTDB_SPLIT_SPECIES {
    label 'process_single'

    container "ghcr.io/labgem/pangbank-wf:merge-split2"

    input:
    path input_genomes
    path taxonomy
    val min_genomes

    output:
    path "gtdb_splits/meta/*.list", optional: true,  emit: genome_list_files
    path "versions.yml", emit: versions

    script:
    """
    find_gtdb_splits.py --input-genomes $input_genomes \
                        --genome-taxonomy $taxonomy \
                        --min-genome-count $min_genomes \
                        --output-directory ./gtdb_splits > find.info

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
