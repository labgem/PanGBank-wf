process FILTER_GENOMES {
    label 'process_single'

    container "ghcr.io/labgem/pangbank-wf:merge-split2"

    input:
    path metadata_file
    path input_genomes
    val genome_min_completeness
    val representative_min_completeness

    output:
    path "genome_quality_filtering.info", emit: genome_quality_filtering_info
    path "input_genomes.filtered.tsv", emit: filtered_genome_input
    path "genome_metadata.filtered.tsv", emit: filtered_genome_metadata
    path "genome_quality_filtering_summary.tsv", emit: genome_quality_filtering_summary
    path "versions.yml", emit: versions

    script:
    """
    filter_genome_on_quality.py --metadata-file $metadata_file \
                        --input-genomes $input_genomes \
                        --genome-min-completeness $genome_min_completeness \
                        --representative-min-completeness $representative_min_completeness \
                        --output-directory ./ > genome_quality_filtering.info

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
