include { METAPANG_INDEX_BANK } from '../../modules/local/metapang_index_bank.nf'
include { METAPANG_INDEX_PANGENOME as METAPANG_INDEX_PANGENOME_SMALL } from '../../modules/local/metapang_index_pangenome.nf'
include { METAPANG_INDEX_PANGENOME as METAPANG_INDEX_PANGENOME_MEDIUM } from '../../modules/local/metapang_index_pangenome.nf'
include { METAPANG_INDEX_PANGENOME as METAPANG_INDEX_PANGENOME_LARGE } from '../../modules/local/metapang_index_pangenome.nf'

workflow METAPANG {
    take:
    pangenomes

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()


    if (params.metapang_build_bank_index) {

        ch_input_bank_index = pangenomes.map { meta, pangenome, genome_file ->
            def paths = genome_file.withInputStream { stream ->
                new java.util.zip.GZIPInputStream(stream).withReader("UTF-8") { reader ->
                    reader.readLines().collect() { line -> line.split ('\t')[1]}.join(' ')
                }
            }
            return "${meta.species} ${paths}\n"
        }.collectFile(name: 'index_bank_input.tsv')


        METAPANG_INDEX_BANK(
            ch_input_bank_index,
            Channel.value(params.metapang_bank_kmer_size),
            Channel.value(params.metapang_bank_scaled)
        )

        ch_bank_index_report = METAPANG_INDEX_BANK.out.stats.collectFile(
            name: 'metapang_index_summary.tsv', keepHeader: true
        )

        ch_versions = METAPANG_INDEX_BANK.out.versions
        ch_multiqc_files = ch_multiqc_files.mix(ch_bank_index_report)
    }

    if (params.metapang_build_pangenome_index) {
        pangenomes_branched = pangenomes.branch { meta, pangenome, genome_file ->
            large: meta.genome_count >= params.large_pangenome_cutoff
            medium: meta.genome_count >= params.large_pangenome_cutoff / 4
            small: true
        }

        METAPANG_INDEX_PANGENOME_LARGE(
            pangenomes_branched.large,
            channel.value(params.metapang_pangenome_kmer_size),
            channel.value(params.metapang_pangenome_annotation)
        )

        METAPANG_INDEX_PANGENOME_MEDIUM(
            pangenomes_branched.medium,
            channel.value(params.metapang_pangenome_kmer_size),
            channel.value(params.metapang_pangenome_annotation)
        )

        METAPANG_INDEX_PANGENOME_SMALL(
            pangenomes_branched.small,
            Channel.value(params.metapang_pangenome_kmer_size),
            Channel.value(params.metapang_pangenome_annotation)
        )

        ch_pangenome_index_report = METAPANG_INDEX_PANGENOME_LARGE.out.stats
            .mix(METAPANG_INDEX_PANGENOME_MEDIUM.out.stats)
            .mix(METAPANG_INDEX_PANGENOME_SMALL.out.stats)
            .collectFile(name: 'metapang_dbg_summary.tsv', keepHeader: true, skip: 1)

        ch_versions = ch_versions.mix(METAPANG_INDEX_PANGENOME_SMALL.out.versions)
        ch_versions = ch_versions.mix(METAPANG_INDEX_PANGENOME_MEDIUM.out.versions)
        ch_versions = ch_versions.mix(METAPANG_INDEX_PANGENOME_LARGE.out.versions)

        ch_multiqc_files = ch_multiqc_files.mix(ch_pangenome_index_report)
    }

    emit:
    multiqc_files = ch_multiqc_files
    versions = ch_versions
}

workflow {
    main:
    ch_versions = channel.empty()

    ch_pangenomes = channel
        .fromPath(params.pangenomes + '/*/input_genomes.tsv.gz')
        .map { f ->
            def n = 0
            new java.util.zip.GZIPInputStream(f.newInputStream()).eachLine { n++ }
            def meta = [
                species     : f.parent.name,
                genome_count: n
            ]
            def h5 = f.resolveSibling('pangenome.h5')
            tuple(meta, h5, f)
        }

    METAPANG(ch_pangenomes)
}
