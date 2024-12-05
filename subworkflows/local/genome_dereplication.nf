/*
 * Dereplicate genomes from species with a lot of genomes
 */

include { SEQFU_STATS_FROM_FILE } from '../../modules/local/seqfu_stats'

include { MASH_SKETCH as MASH_SKETCH_GENOMES } from '../../modules/nf-core/mash/sketch/main'
include { SORT_GENOMES } from '../../modules/local/sort_genomes'
include { MASH_DIST_TO_PHYLIP } from '../../modules/local/mash_distance_in_phylip'

include { QUICKTREE } from '../../modules/local/quicktree'
include { GENOME_SELECTION_FROM_TREE } from '../../modules/local/genome_selection_from_tree'
include { ANY2FASTA } from '../../modules/local/any2fasta'

include { FORMAT_INPUT_GENOMES } from '../../modules/local/format_input_genomes'





workflow GENOME_DEREPLICATION {
    take:
    ch_species_to_dereplicate   // [ [ meta] , fasta    ],

    main:
    ch_versions = Channel.empty()
    ch_genome_count_cutoff = Channel.value(params.dereplication_genome_cutoff)

    // ch_species_to_dereplicate.view()

    ch_species_branched = ch_species_to_dereplicate.branch{meta, genome_file ->
                                                            annotation_input : meta.file_type == "annotation"
                                                            fasta_input : true}
    // ch_species_branched.annotation_input.view { v -> "annotation_input: ${v}" }
    // ch_species_branched.fasta_input.view { v -> "fasta_input: ${v}" }

    ANY2FASTA(ch_species_branched.annotation_input)

    // ANY2FASTA.out.genome_path_fasta.view()

    ch_species_to_fasta_input_files = ch_species_branched.fasta_input.concat(ANY2FASTA.out.genome_path_fasta)

    // mash can take a file with genome paths in it
    // need to convert of genome_name path file to a single path list file
    ch_species_to_path_file = ch_species_to_fasta_input_files.splitCsv(elem:1, sep:"\t", header: ["name", "path"])
                                                            .collectFile(newLine: true){ meta, genome -> ["${meta.species}", genome.path]}
                                                            .map{path_file -> [["id":path_file.name], path_file]}



    SEQFU_STATS_FROM_FILE(ch_species_to_path_file)

    SORT_GENOMES(SEQFU_STATS_FROM_FILE.out.stats)

    MASH_SKETCH_GENOMES(ch_species_to_path_file)

    ch_species_sketch_genome_list = MASH_SKETCH_GENOMES.out.mash.concat(SORT_GENOMES.out.sorted_genomes_list)
                                                            .groupTuple()
                                                            .map {meta, files ->
                                                            def mash_sketch = files[0]
                                                            def genome_list = files[1]
                                                            [meta, mash_sketch, genome_list]}


    MASH_DIST_TO_PHYLIP(ch_species_sketch_genome_list)

    QUICKTREE(MASH_DIST_TO_PHYLIP.out.phylip_matrix)

    ch_sp_and_genome_name_to_path =  ch_species_to_fasta_input_files.map{ meta, genome_file -> [["id":meta.species], genome_file] }
    ch_sp_and_fasta_to_orginal_path = ANY2FASTA.out.fasta_to_orginal_path.map{meta, file -> [["id":meta.species], file]}


    ch_sp_tree_and_sorted_genomes = QUICKTREE.out.tree.concat(SORT_GENOMES.out.sorted_genomes_list)
                                            .groupTuple()
                                            .map {meta, files ->
                                                    def tree = files[0]
                                                    def sorted_genomes = files[1]
                                                    [meta, tree, sorted_genomes]}



    // ch_sp_tree_genome_list_and_name_to_path = QUICKTREE.out.tree.concat(SORT_GENOMES.out.sorted_genomes_list, ch_sp_and_genome_name_to_path, ch_sp_and_original_path_to_fasta)
    //                                         .groupTuple()
    //                                         .map {meta, files ->
    //                                                 def tree = files[0]
    //                                                 def sorted_genomes = files[1]
    //                                                 def genome_name_to_path = files[2]
    //                                                 def fasta_to_original_input = files.size() <= 3 ? file("NO_FILE") : files[3]
    //                                                 [meta, tree, sorted_genomes, genome_name_to_path, fasta_to_original_input]}


    GENOME_SELECTION_FROM_TREE(ch_sp_tree_and_sorted_genomes, ch_genome_count_cutoff)

    GENOME_SELECTION_FROM_TREE.out.selected_genomes.view{i -> "AAAAAAAA $i "}

    ch_sp_selected_genome_to_name_and_original_path = GENOME_SELECTION_FROM_TREE.out.selected_genomes.concat(ch_sp_and_genome_name_to_path, ch_sp_and_fasta_to_orginal_path)
                                                    .groupTuple()
                                                    .map {meta, files ->
                                                            def selected_genomes = files[0]
                                                            def genome_name_to_path = files[1]
                                                            def fasta_to_original_input = files.size() <= 2 ? file("NO_FILE") : files[2]
                                                            [meta, selected_genomes, genome_name_to_path, fasta_to_original_input]}


    ch_sp_selected_genome_to_name_and_original_path.view{i -> "BBBBBBBB $i "}

    FORMAT_INPUT_GENOMES(ch_sp_selected_genome_to_name_and_original_path)

    // ch_cluster_compo_and_phylip_matrix = GENOME_SELECTION_FROM_TREE.out.cluster_composition.concat(MASH_DIST_TO_PHYLIP.out.phylip_matrix)
    //                                                             .groupTuple(size:2)
    //                                                             .map {meta, cluster_compo_and_phylip_matrix
    //                                                                     -> [meta, cluster_compo_and_phylip_matrix[0], cluster_compo_and_phylip_matrix[1]]}
    //                                                             // .view  { v -> "ch_cluster_compo_and_phylip_matrix: ${v}" }
    // // CLUSTER_STAT()

    // // Build final channel to emit
    ch_sp_to_original_meta = ch_species_to_dereplicate.map { meta, genome_file -> [["id":meta.species], meta]}
    ch_meta_and_selected_genomes = ch_sp_to_original_meta.concat(FORMAT_INPUT_GENOMES.out.genome_selection_ppanggo_input)
                                                    .groupTuple(size:2)
                                                    .map {sp, meta_and_selected_genomes ->
                                                        def meta = meta_and_selected_genomes[0]
                                                        def selected_genomes =  meta_and_selected_genomes[1]
                                                        meta.genomes_count = selected_genomes.countLines()
                                                        [meta, selected_genomes]}
                                                    .view { i -> "CCCCCC $i"}


    // ch_versions.mix( MASH_SKETCH_GENOMES.out.versions )

    emit:
    dereplicated_genomes = ch_meta_and_selected_genomes
    versions = ch_versions

}
