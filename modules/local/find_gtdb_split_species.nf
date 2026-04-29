process FIND_GTDB_SPLIT_SPECIES {
    label 'process_single'

    container "ghcr.io/labgem/pangbank-wf:merge-split2"

    input:
    path metadata_file
    val genome_min_checkm
    val genome_min_checkm2
    val representative_min_checkm
    val representative_min_checkm2
    val min_genomes

    output:
    path "gtdb_splits/meta/*.list", emit: genome_list_files
    path "versions.yml", emit: versions

    script:
    """
    find_gtdb_splits.py --metadata-file $metadata_file \
                        --genome-min-checkm $genome_min_checkm \
                        --genome-min-checkm2 $genome_min_checkm2 \
                        --representative-min-checkm $representative_min_checkm \
                        --representative-min-checkm2 $representative_min_checkm2 \
                        --min-genome-count $min_genomes \
                        --output-directory ./gtdb_splits

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
