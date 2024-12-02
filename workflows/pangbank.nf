/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_ppanggolin_config          = Channel.fromPath("$projectDir/assets/ppanggolin_config.yml", checkIfExists: true)
// ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
// ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
// ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
// ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { GENOME_DEREPLICATION } from '../subworkflows/local/genome_dereplication'

//
// MODULE: Local modules
//
include { PARSE_GENOMES_AND_TAXONOMY                      } from '../modules/local/parse_genomes_and_taxonomy'
include { PPANGGOLIN_ALL                                  } from '../modules/local/ppanggolin/all'
include { PPANGGOLIN_FASTA                                } from '../modules/local/ppanggolin/fasta'
include { GATHER_PANGENOME_INFO                           } from '../modules/local/gather_pangenome_infos'
include { MD5SUM_ON_FILES                                 } from '../modules/local/md5sum_on_list_of_files'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// MODULE: Installed directly from nf-core/modules

include { MASH_SKETCH                                     } from '../modules/nf-core/mash/sketch/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PANGBANK {

    main:
    ch_versions = Channel.empty()

    ch_min_genomes = Channel.value(params.min_genomes)

    ch_input_genomes = manage_input_genomes(file(params.genomes))

    // PREPARE SPECIES: Check species that have enough genome to build a pangenome
    PARSE_GENOMES_AND_TAXONOMY (
        ch_input_genomes,
        file(params.taxonomy),
        ch_min_genomes
    )

    ch_versions = ch_versions.mix(PARSE_GENOMES_AND_TAXONOMY.out.versions)

    ch_ppanggo_inputs_meta = PARSE_GENOMES_AND_TAXONOMY.out.ppanggo_inputs.flatten()
                                                .map { create_ppanggo_input_channel(it) }

    // ch_ppanggo_inputs_meta.view()
    // ch_species_to_dereplicate = ch_ppanggo_inputs_meta.filter{meta, genome_file -> meta.genomes_count > params.dereplication_genome_cutoff}

    // ch_ppanggo_inputs_meta.view()
    ch_species_branched = ch_ppanggo_inputs_meta.branch{meta, genome_file ->
                                                            to_dereplicate : meta.genomes_count > params.dereplication_genome_cutoff
                                                            other : true}

    // ch_species_branched.view()
    // ch_species_to_dereplicate.view()
    // ch_species_branched.to_dereplicate.view()
    // dereplication_genome_cutoff
    GENOME_DEREPLICATION(ch_species_branched.to_dereplicate)

    PPANGGOLIN_ALL(ch_ppanggo_inputs_meta, ch_ppanggolin_config.toList())

    PPANGGOLIN_FASTA(PPANGGOLIN_ALL.out.pangenome)

    ch_fasta_list_file = PPANGGOLIN_FASTA.out.persistent_families_fasta.map{meta, fasta -> fasta.path}
                                                                    .collectFile(name: 'persistent_fasta_list.txt', newLine: true)
                                                                    .map{file -> [[id:"families_persistent_all.msh"], file]}


    MASH_SKETCH(ch_fasta_list_file)



    MD5SUM_ON_FILES(ch_ppanggo_inputs_meta)

    ch_pangenome_infos = PPANGGOLIN_ALL.out.pangenome_info.collect()

    GATHER_PANGENOME_INFO(ch_pangenome_infos)

}

def manage_input_genomes(input_genomes_file) {

    // This function process input genomes
    // If paths in the file are URLs, it downloads the files and stores them in a temporary location,
    // generating a new genome list file mapping genome names to local paths.

    log.info "Processing input genomes file: ${input_genomes_file}"

    // Read the first line of the file to get the genome file path or URL
    def first_line = ""
    input_genomes_file.eachLine { line ->
        first_line = line
        return false  // Stop after reading the first line
    }
    // This is cleaner bu produced an error in github CI with latest nextflow version
    // def first_line = input_genomes_file.withReader('UTF-8') { it.readLine() }
    // log.info "First line from genome file: ${first_line}"

    // Extract the genome file path or URL (assuming second column after splitting by tab)
    def genome_file_str = first_line.split('\t')[1]

    // List of supported URL schemes
    def url_schemes = ["https", "http", "ftp"]

    // Regular expression pattern to extract the file extension (excluding 'gz')
    def extension_pattern = ~/.*((\.[a-yA-Y]+)(\.gz)?)$/

    // Check if the genome file path is a URL
    if (url_schemes.contains(file(genome_file_str).getScheme())) {
        // If it's a URL, download the files and map them to local paths
        log.info "Downloading genome files from URLs"

        def genome_file_channel = channel.fromPath(input_genomes_file)
            .splitCsv(header: ['name', 'path'], sep: "\t")
            .map { row ->
                def file_extension = (row.path =~ extension_pattern)[0][1]
                ["file_name": "${row.name}${file_extension}", "path": row.path]
            }
            .collectFile { row ->
                ["${row.file_name}", file(row.path).bytes]
            }
            .map { path ->
                def file_extension = (path.name =~ extension_pattern)[0][1]
                ["name": "${path.name - file_extension}", "path": path]
            }
            .collectFile(name: 'input_genomes.txt', newLine: true) { row ->
                ['input_genomes.txt', "${row.name}\t${row.path}"]
            }

        return genome_file_channel

    } else {
        // If it's a local file, process the file paths directly
        log.info "Processing local genome file paths"

        def genome_file_channel = channel.fromPath(input_genomes_file)
            .splitCsv(header: ['name', 'path'], sep: "\t")
            .collectFile(name: 'input_genomes.txt', newLine: true) { row ->
                ['input_genomes.txt', "${row.name}\t${file(row.path)}"]
            }

        return genome_file_channel
    }
}



def create_ppanggo_input_channel(input_file) {

    // create meta map
    def meta = [:]


    meta.species = input_file.getSimpleName()
    meta.genomes_count = input_file.countLines()

    // Getting the extension of the first genome file
    def first_line = input_file.withReader { it.readLine() }
    def first_genome_file = first_line.split('\t')[1]

    def extension_patern = ~/.*(\.[a-yA-Y]+)(\.gz)?$/ // A-Y to exclude Z to not cacth gz if exists.
    def genome_extension = (first_genome_file =~ extension_patern)[0][1].toLowerCase()
    def annotation_extensions = params.annotation_extensions.split(';')
    def fasta_extensions = params.fasta_extensions.split(';')

    if (annotation_extensions.contains(genome_extension)){
        meta.file_type = "annotation"
    } else if (fasta_extensions.contains(genome_extension) ) {
        meta.file_type = "fasta"
    }
    else {
        exit 1, """
        ERROR: Please check input genomes -> Genome file (${first_genome_file}) does have an unexpected extension: $genome_extension
        Possible value for annotation files: $annotation_extensions\
        Fasta files: $fasta_extensions
        """
    }
    def input_meta = [ meta, input_file]

    return input_meta
}
