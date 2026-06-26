include { FIND_GTDB_SPLIT_SPECIES } from '../../modules/local/find_gtdb_split_species.nf'
include { MERGE_GTDB_SPLIT_SPECIES as MERGE_GTDB_SPLIT_SPECIES_SMALL } from '../../modules/local/merge_gtdb_split_species.nf'
include { MERGE_GTDB_SPLIT_SPECIES as MERGE_GTDB_SPLIT_SPECIES_MEDIUM } from '../../modules/local/merge_gtdb_split_species.nf'
include { MERGE_GTDB_SPLIT_SPECIES as MERGE_GTDB_SPLIT_SPECIES_LARGE } from '../../modules/local/merge_gtdb_split_species.nf'

workflow GTDB_SPLIT_SPECIES {
    take:
    input_genomes
    genome_fasta
    taxonomy
    min_genomes


    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    FIND_GTDB_SPLIT_SPECIES(
        input_genomes,
        taxonomy,
        min_genomes
    )
    ch_versions = ch_versions.mix(FIND_GTDB_SPLIT_SPECIES.out.versions)

    ch_split_species = FIND_GTDB_SPLIT_SPECIES.out.genome_list_files
                                                .flatten()
                                                .map { genome_list_file -> tuple([id: genome_list_file.baseName,
                                                                                  genomes_count: genome_list_file.countLines()],
                                                                                  genome_list_file) }




    ch_split_species_branched = ch_split_species.branch { meta, _genome_file ->
        large: meta.genomes_count >= params.large_pangenome_cutoff * 4
        medium: meta.genomes_count >= params.large_pangenome_cutoff / 5
        small: true
    }


    MERGE_GTDB_SPLIT_SPECIES_LARGE(ch_split_species_branched.large,
                                    genome_fasta,
                                    params.gtdb_merge_ani_threshold,
                                    params.gtdb_merge_af_threshold,
                                    params.gtdb_merge_species_ani_stat)
    MERGE_GTDB_SPLIT_SPECIES_MEDIUM(ch_split_species_branched.medium,
                                    genome_fasta,
                                    params.gtdb_merge_ani_threshold,
                                    params.gtdb_merge_af_threshold,
                                    params.gtdb_merge_species_ani_stat)
    MERGE_GTDB_SPLIT_SPECIES_SMALL(ch_split_species_branched.small,
                                    genome_fasta,
                                    params.gtdb_merge_ani_threshold,
                                    params.gtdb_merge_af_threshold,
                                    params.gtdb_merge_species_ani_stat)

    ch_versions = ch_versions.mix(MERGE_GTDB_SPLIT_SPECIES_SMALL.out.versions)
    ch_versions = ch_versions.mix(MERGE_GTDB_SPLIT_SPECIES_MEDIUM.out.versions)
    ch_versions = ch_versions.mix(MERGE_GTDB_SPLIT_SPECIES_LARGE.out.versions)

    split_clusters = MERGE_GTDB_SPLIT_SPECIES_SMALL.out.split_clusters
                        .concat(MERGE_GTDB_SPLIT_SPECIES_MEDIUM.out.split_clusters,
                                MERGE_GTDB_SPLIT_SPECIES_LARGE.out.split_clusters)

    genome_clusters = MERGE_GTDB_SPLIT_SPECIES_SMALL.out.genome_clusters
                        .concat(MERGE_GTDB_SPLIT_SPECIES_MEDIUM.out.genome_clusters,
                                MERGE_GTDB_SPLIT_SPECIES_LARGE.out.genome_clusters)

    species_pair_summary = MERGE_GTDB_SPLIT_SPECIES_SMALL.out.species_pair_summary
                        .concat(MERGE_GTDB_SPLIT_SPECIES_MEDIUM.out.species_pair_summary,
                                MERGE_GTDB_SPLIT_SPECIES_LARGE.out.species_pair_summary)


    ch_split_clusters = split_clusters.collectFile(name: 'split_clusters.tsv', storeDir: "${params.outdir}/genome_preprocessing/merge_gtdb_split/").ifEmpty(file("$projectDir/assets/NO_FILE_3"))
    ch_genome_clusters = genome_clusters.collectFile(name: 'genome_clusters.tsv', storeDir: "${params.outdir}/genome_preprocessing/merge_gtdb_split/")

    _ch_species_pair_summary = species_pair_summary.collectFile(name: 'species_pair_summary.tsv', storeDir: "${params.outdir}/genome_preprocessing/merge_gtdb_split/", keepHeader: true)

    ch_gtdb_merge_summary = MERGE_GTDB_SPLIT_SPECIES_SMALL.out.merge_summary
                            .concat(MERGE_GTDB_SPLIT_SPECIES_MEDIUM.out.merge_summary,
                                    MERGE_GTDB_SPLIT_SPECIES_LARGE.out.merge_summary)
                            .collectFile(skip: 1, keepHeader: true, name: 'gtdb_merge_summary.tsv')

    ch_multiqc_files = ch_multiqc_files.mix(ch_gtdb_merge_summary)

    emit:
    split_clusters = ch_split_clusters
    genome_clusters = ch_genome_clusters
    multiqc_files = ch_multiqc_files
    versions = ch_versions

}
