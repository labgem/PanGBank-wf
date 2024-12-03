/*
 * Dereplicate genomes from species with a lot of genomes
 */

include { SEQFU_STATS_FROM_FILE } from '../../modules/local/seqfu_stats'

include { MASH_SKETCH as MASH_SKETCH_GENOMES } from '../../modules/nf-core/mash/sketch/main'
include { SORT_GENOMES } from '../../modules/local/sort_genomes'
include { MASH_DIST_TO_PHYLIP } from '../../modules/local/mash_distance_in_phylip'

include { QUICKTREE } from '../../modules/local/quicktree'


workflow GENOME_DEREPLICATION {
    take:
    ch_species_to_dereplicate   // [ [ meta] , fasta    ],

    main:
    ch_versions = Channel.empty()

    // mash can take a file with genome paths in it
    // need to convert of genome_name path file to a single path list file
    ch_species_to_path_file = ch_species_to_dereplicate.splitCsv(elem:1, sep:"\t", header: ["name", "path"])
                            .collectFile(newLine: true){ meta, genome -> ["${meta.species}", genome.path]}
                            .map{path_file -> [["id":path_file.name], path_file]}


    // SEQFU_STATS from nf core
    // ch_meta_to_single_genome = ch_species_to_dereplicate.splitCsv(elem:1, sep:"\t", header: ["name", "path"])
    //                             .map{meta, content -> [meta, content.path]}

    // SEQFU_STATS(ch_meta_to_single_genome)

    // ch_meta_genome_stat = SEQFU_STATS.out.stats
    //                         .collectFile(keepHeader:true){ meta, genome_stat -> ["${meta.species}", genome_stat]}
    //                         .map{path_file -> [["id":path_file.name], path_file]}

    // ch_meta_genome_stat.view()

    // SORT_GENOMES(ch_meta_genome_stat)
    //SEQFU_STATS.out.stats.view()

    // input could be splitted to run faster
    SEQFU_STATS_FROM_FILE(ch_species_to_path_file)

    SORT_GENOMES(SEQFU_STATS_FROM_FILE.out.stats)

    MASH_SKETCH_GENOMES(ch_species_to_path_file)
    species_sketch_genome_list = MASH_SKETCH_GENOMES.out.mash.concat(SORT_GENOMES.out.sorted_genomes_list)
                                                            .view  { v -> "concat: ${v}" }
                                                            .groupTuple()
                                                            .map {meta, sketch_and_genome_list -> [meta, sketch_and_genome_list[0], sketch_and_genome_list[1]]}
                                                            .view  { v -> "groupTuple: ${v}" }


    MASH_DIST_TO_PHYLIP(species_sketch_genome_list)

    QUICKTREE(MASH_DIST_TO_PHYLIP.out.phylip_matrix)


    // GENOME_SELECTION_FROM_TREE(MASH_DIST_TO_PHYLIP.out.phylip_matrix)

    // ch_versions.mix( MASH_SKETCH_GENOMES.out.versions )

    emit:
    versions = ch_versions

}
