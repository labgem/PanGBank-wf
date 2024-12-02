/*
 * Dereplicate genomes from species with a lot of genomes
 */

// include { GENOMAD_DOWNLOAD   } from '../../modules/nf-core/genomad/download/main'
// include { GENOMAD_ENDTOEND   } from '../../modules/nf-core/genomad/endtoend/main'
include { SEQFU_STATS } from '../../modules/nf-core/seqfu/stats/main'
include { MASH_SKETCH as MASH_SKETCH_GENOMES } from '../../modules/nf-core/mash/sketch/main'
include { SORT_GENOMES } from '../../modules/local/sort_genomes'



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

    ch_meta_to_single_genome = ch_species_to_dereplicate.splitCsv(elem:1, sep:"\t", header: ["name", "path"])
                                .map{meta, content -> [meta, content.path]}

    SEQFU_STATS(ch_meta_to_single_genome)

    ch_meta_genome_stat = SEQFU_STATS.out.stats
                            .collectFile(keepHeader:true){ meta, genome_stat -> ["${meta.species}", genome_stat]}
                            .map{path_file -> [["id":path_file.name], path_file]}

    ch_meta_genome_stat.view()

    SORT_GENOMES(ch_meta_genome_stat)
    //SEQFU_STATS.out.stats.view()

    MASH_SKETCH_GENOMES(ch_species_to_path_file)

    // MASH_DIST_TO_PHYLIP(MASH_SKETCH_GENOMES.out.mash)

    ch_versions.mix( MASH_SKETCH_GENOMES.out.versions )

    emit:
    versions = ch_versions

}
