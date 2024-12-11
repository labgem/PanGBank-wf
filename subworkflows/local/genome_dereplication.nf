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
include { CLUSTER_STAT } from '../../modules/local/cluster_stat.nf'
include { CLUSTER_PLOT } from '../../modules/local/cluster_plot.nf'

include { mergeText } from 'plugin/nf-boost'

workflow GENOME_DEREPLICATION {
    take:
    ch_species_to_dereplicate // Channel with species metadata and genome path files.

    main:
    ch_versions = Channel.empty()
    ch_genome_count_cutoff = Channel.value(params.dereplication_threshold)

    if (params.reference_genomes) {
        ch_reference_genomes = file(params.reference_genomes, checkIfExists: true)
    }
    else {
        ch_reference_genomes = []
    }

    // Branch the input channel based on the type of input file (annotation or fasta).
    ch_species_branched = ch_species_to_dereplicate.branch { meta, genome_file ->
        annotation_input: meta.file_type == "annotation"
        fasta_input: true
    }

    ch_sp_annotation_input_split = ch_species_branched.annotation_input.splitText(elem: 1, by: 500, file: true)


    // Convert annotation files to FASTA format if required.
    ANY2FASTA(ch_sp_annotation_input_split)
    ch_versions = ch_versions.mix(ANY2FASTA.out.versions)

    ch_sp_and_fasta_to_original_path = ANY2FASTA.out.fasta_to_orginal_path
        .collectFile(newLine: false) { meta, content -> ["${meta.species}.fasta_to_orginal_path", content] }
        .map { path_file -> [["id": path_file.baseName], path_file] }

    ch_sp_and_genome_path_fasta = ANY2FASTA.out.genome_path_fasta
        .collectFile(newLine: false) { meta, content -> ["${meta.species}.fasta_input_file", content] }
        .map { path_file -> [["id": path_file.baseName], path_file] }


    // Combine fasta input files from both branches (direct and converted).
    ch_species_to_fasta_input_files = ch_species_branched.fasta_input
        .map { meta, file -> [["id": meta.species], file] }
        .concat(ch_sp_and_genome_path_fasta)

    // Generate a single path list file for all genomes within a species.
    ch_species_to_path_file = ch_species_to_fasta_input_files
        .splitCsv(elem: 1, sep: "\t", header: ["name", "path"])
        .collectFile(newLine: true) { meta, genome -> ["${meta.id}.fasta_input_file", genome.path] }
        .map { path_file -> [["id": path_file.baseName], path_file] }



    // Calculate genome sequence metrics using seqfu.
    SEQFU_STATS_FROM_FILE(ch_species_to_path_file)
    ch_versions = ch_versions.mix(SEQFU_STATS_FROM_FILE.out.versions)

    // Sort genomes based on calculated metrics for quality prioritization.
    SORT_GENOMES(SEQFU_STATS_FROM_FILE.out.stats)
    ch_versions = ch_versions.mix(SORT_GENOMES.out.versions)

    // Compute pairwise Mash distances for all genomes.
    MASH_SKETCH_GENOMES(ch_species_to_path_file)
    ch_versions = ch_versions.mix(MASH_SKETCH_GENOMES.out.versions)

    // Combine Mash sketches with sorted genome lists for further processing.
    ch_species_sketch_genome_list = MASH_SKETCH_GENOMES.out.mash
        .concat(SORT_GENOMES.out.sorted_genomes_list)
        .groupTuple()
        .map { meta, files ->
            def mash_sketch = files[0]
            def genome_list = files[1]
            [meta, mash_sketch, genome_list]
        }


    // Convert Mash distances to a PHYLIP matrix format.
    MASH_DIST_TO_PHYLIP(ch_species_sketch_genome_list)
    ch_versions = ch_versions.mix(MASH_DIST_TO_PHYLIP.out.versions)

    // Build a phylogenetic tree from the distance matrix.
    QUICKTREE(MASH_DIST_TO_PHYLIP.out.phylip_matrix)
    ch_versions = ch_versions.mix(QUICKTREE.out.versions)

    // Combine the tree with sorted genomes for clustering.
    ch_sp_tree_and_sorted_genomes = QUICKTREE.out.tree
        .concat(SORT_GENOMES.out.sorted_genomes_list)
        .groupTuple()
        .map { meta, files ->
            def tree = files[0]
            def sorted_genomes = files[1]
            [meta, tree, sorted_genomes]
        }

    // Select genome clusters from the phylogenetic tree and prioritize based on quality.
    GENOME_SELECTION_FROM_TREE(ch_sp_tree_and_sorted_genomes, ch_genome_count_cutoff)
    ch_versions = ch_versions.mix(GENOME_SELECTION_FROM_TREE.out.versions)

    // Prepare input files for PPanGGOLiN by mapping selected genome paths to their original inputs.
    ch_sp_selected_genome_to_name_and_original_path = GENOME_SELECTION_FROM_TREE.out.selected_genomes
        .concat(ch_species_to_fasta_input_files, ch_sp_and_fasta_to_original_path)
        .groupTuple()
        .map { meta, files ->
            def selected_genomes = files[0]
            def genome_name_to_path = files[1]
            def fasta_to_original_input = files.size() <= 2 ? file("NO_FILE") : files[2]
            [meta, selected_genomes, genome_name_to_path, fasta_to_original_input]
        }
    FORMAT_INPUT_GENOMES(ch_sp_selected_genome_to_name_and_original_path, ch_reference_genomes)
    ch_versions = ch_versions.mix(FORMAT_INPUT_GENOMES.out.versions)


    // Compute cluster statistics and generate a cluster composition file.
    ch_cluster_compo_and_phylip_matrix = GENOME_SELECTION_FROM_TREE.out.cluster_composition
        .concat(MASH_DIST_TO_PHYLIP.out.phylip_matrix)
        .groupTuple(size: 2)
        .map { meta, cluster_compo_and_phylip_matrix ->
            [meta, cluster_compo_and_phylip_matrix[0], cluster_compo_and_phylip_matrix[1]]
        }
    CLUSTER_STAT(ch_cluster_compo_and_phylip_matrix)
    ch_versions = ch_versions.mix(CLUSTER_STAT.out.versions)

    // Generate plots from cluster statistics.
    ch_cluster_stat = CLUSTER_STAT.out.cluster_stat
        .map { meta, cluster_stat -> cluster_stat }
        .collectFile(skip: 1, keepHeader: true, name: 'cluster_stat.tsv')
    ch_distance_count = CLUSTER_STAT.out.distance_count
        .map { meta, distance_count -> distance_count }
        .collectFile(skip: 1, keepHeader: true, name: 'distance_count.tsv')

    CLUSTER_PLOT(ch_cluster_stat, ch_distance_count)
    ch_versions = ch_versions.mix(CLUSTER_PLOT.out.versions)


    // Merge updated metadata with the selected genomes and finalize output.
    ch_sp_to_original_meta = ch_species_to_dereplicate.map { meta, genome_file -> [["id": meta.species], meta] }
    ch_meta_and_selected_genomes = ch_sp_to_original_meta
        .concat(FORMAT_INPUT_GENOMES.out.genome_selection_ppanggo_input)
        .groupTuple(size: 2)
        .map { sp, meta_and_selected_genomes ->
            def meta = meta_and_selected_genomes[0]
            def selected_genomes = meta_and_selected_genomes[1]
            meta.genomes_count = selected_genomes.countLines()
            [meta, selected_genomes]
        }

    emit:
    dereplicated_genomes = ch_meta_and_selected_genomes
    versions = ch_versions
}
