process FIND_GTDB_SPLIT_SPECIES {
    label 'process_single'

    container "ghcr.io/labgem/pangbank-wf:merge-split2"

    input:
    path metadata_file
    path input_genomes
    val genome_min_completeness
    val representative_min_completeness
    val min_genomes

    output:
    path "gtdb_splits/meta/*.list", emit: genome_list_files
    path "versions.yml", emit: versions

    script:
    """
    find_gtdb_splits.py --metadata-file $metadata_file \
                        --used-genomes $input_genomes \
                        --genome-min-completeness $genome_min_completeness \
                        --representative-min-completeness $representative_min_completeness \
                        --min-genome-count $min_genomes \
                        --output-directory ./gtdb_splits > find.info

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
