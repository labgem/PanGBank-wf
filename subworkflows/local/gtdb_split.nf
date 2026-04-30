include { FIND_GTDB_SPLIT_SPECIES } from '../../modules/local/find_gtdb_split_species.nf'
include { MERGE_GTDB_SPLIT_SPECIES } from '../../modules/local/merge_gtdb_split_species.nf'

workflow GTDB_SPLIT_SPECIES {
    take:
    input_genomes
    genome_fasta
    taxonomy
    min_genomes


    main:
    ch_versions = channel.empty()

    FIND_GTDB_SPLIT_SPECIES(
        input_genomes,
        taxonomy,
        min_genomes
    )
    ch_versions = ch_versions.mix(FIND_GTDB_SPLIT_SPECIES.out.versions)

    ch_split_species = FIND_GTDB_SPLIT_SPECIES.out.genome_list_files

    MERGE_GTDB_SPLIT_SPECIES(
        ch_split_species.flatten().map { f -> tuple([id: f.baseName], f) },
        genome_fasta.collect(),
        params.gtdb_merge_ani_threshold
    )
    ch_versions = ch_versions.mix(MERGE_GTDB_SPLIT_SPECIES.out.versions)

    split_clusters = MERGE_GTDB_SPLIT_SPECIES.out.split_clusters.collectFile(name: 'split_clusters.tsv').ifEmpty(file("$projectDir/assets/NO_FILE_3"))

    emit:
    split_clusters = split_clusters
    genome_clusters = MERGE_GTDB_SPLIT_SPECIES.out.genome_clusters
    versions = ch_versions

}
