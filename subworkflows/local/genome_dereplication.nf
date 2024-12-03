/*
 * Dereplicate genomes from species with a lot of genomes
 */

include { SEQFU_STATS_FROM_FILE } from '../../modules/local/seqfu_stats'

include { MASH_SKETCH as MASH_SKETCH_GENOMES } from '../../modules/nf-core/mash/sketch/main'
include { SORT_GENOMES } from '../../modules/local/sort_genomes'
include { MASH_DIST_TO_PHYLIP } from '../../modules/local/mash_distance_in_phylip'

include { QUICKTREE } from '../../modules/local/quicktree'
include { GENOME_SELECTION_FROM_TREE } from '../../modules/local/genome_selection_from_tree'

workflow GENOME_DEREPLICATION {
    take:
    ch_species_to_dereplicate   // [ [ meta] , fasta    ],

    main:
    ch_versions = Channel.empty()
    ch_genome_count_cutoff = Channel.value(params.dereplication_genome_cutoff)


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
                                                            .groupTuple()
                                                            .map {meta, sketch_and_genome_list -> [meta, sketch_and_genome_list[0], sketch_and_genome_list[1]]}


    MASH_DIST_TO_PHYLIP(species_sketch_genome_list)

    QUICKTREE(MASH_DIST_TO_PHYLIP.out.phylip_matrix)

    ch_sp_and_genome_name_to_path =  ch_species_to_dereplicate.map{ meta, genome_file -> [["id":meta.species], genome_file] }

    ch_sp_tree_genome_list_and_name_to_path = QUICKTREE.out.tree.concat(SORT_GENOMES.out.sorted_genomes_list, ch_sp_and_genome_name_to_path)
                                            .groupTuple(size:3)
                                            .view  { v -> "-->groupTuple: ${v}" }
                                            .map {meta, tree_list_and_name_to_path
                                                    -> [meta, tree_list_and_name_to_path[0], tree_list_and_name_to_path[1], tree_list_and_name_to_path[2]]}
                                            .view  { v -> "map: ${v}" }


    GENOME_SELECTION_FROM_TREE(ch_sp_tree_genome_list_and_name_to_path, ch_genome_count_cutoff)

    // ch_versions.mix( MASH_SKETCH_GENOMES.out.versions )

    emit:
    versions = ch_versions

}
