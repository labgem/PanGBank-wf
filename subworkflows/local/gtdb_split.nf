include { FIND_GTDB_SPLIT_SPECIES } from '../../modules/local/find_gtdb_split_species.nf'
include { MERGE_GTDB_SPLIT_SPECIES } from '../../modules/local/merge_gtdb_split_species.nf'

workflow GTDB_SPLIT_SPECIES {
    take:
    input_genomes
    genome_fasta

    main:
    ch_versions = channel.empty()

    FIND_GTDB_SPLIT_SPECIES(
        file(params.genome_metadata),
        input_genomes,
        params.genome_min_checkm,
        params.representative_min_checkm,
        params.min_genomes
    )
    ch_versions = ch_versions.mix(FIND_GTDB_SPLIT_SPECIES.out.versions)

    ch_split_species = FIND_GTDB_SPLIT_SPECIES.out.genome_list_files

    MERGE_GTDB_SPLIT_SPECIES(
        ch_split_species.flatten().map { f -> tuple([id: f.baseName], f) },
        genome_fasta.collect(),
        params.gtdb_merge_threshold
    )
    ch_versions = ch_versions.mix(MERGE_GTDB_SPLIT_SPECIES.out.versions)

    emit:
    split_clusters = MERGE_GTDB_SPLIT_SPECIES.out.split_clusters
    genome_clusters = MERGE_GTDB_SPLIT_SPECIES.out.genome_clusters
    versions = ch_versions

}
