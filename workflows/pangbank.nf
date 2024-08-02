/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowPangbank.initialise(params, log)

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
// include { INPUT_CHECK } from '../subworkflows/local/input_check'

//
// MODULE: Local modules
//
include { PARSE_GENOMES_AND_TAXONOMY                      } from '../modules/local/parse_genomes_and_taxonomy'
include { PPANGGOLIN                                      } from '../modules/local/ppanggolin'
include { GATHER_PANGENOME_INFO                           } from '../modules/local/gather_pangenome_infos'
include { MD5SUM_ON_FILES                                 } from '../modules/local/md5sum_on_list_of_files'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
// include { FASTQC                      } from '../modules/nf-core/fastqc/main'
// include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
// include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
// def multiqc_report = []

workflow PANGBANK {

    ch_versions = Channel.empty()

    ch_min_genomes = Channel.value(params.min_genomes)

    // PREPARE SPECIES: Check species that have enough genome to build a pangenome
    PARSE_GENOMES_AND_TAXONOMY (
        file(params.genomes),
        file(params.taxonomy),
        ch_min_genomes
    )

    ch_versions = ch_versions.mix(PARSE_GENOMES_AND_TAXONOMY.out.versions)

    ch_ppanggo_inputs_meta = PARSE_GENOMES_AND_TAXONOMY.out.ppanggo_inputs.flatten()
                                                .map { create_ppanggo_input_channel(it) }

    ch_ppanggo_inputs_meta.view()

    PPANGGOLIN(ch_ppanggo_inputs_meta, ch_ppanggolin_config.toList())

    MD5SUM_ON_FILES(ch_ppanggo_inputs_meta)


    ch_pangenome_infos = PPANGGOLIN.out.pangenome_info.collect()

    GATHER_PANGENOME_INFO(ch_pangenome_infos)

}

def create_ppanggo_input_channel(input_file) {

    // def annotation_exts = [".gbff", ".gff", ".gb"];
    // def fasta_exts = [".fna", ".fasta", ".fa"];
    // create meta map
    def meta = [:]
    def genome_files = []

    meta.species = input_file.getSimpleName()
    meta.genomes_count = input_file.countLines()

    // need to extract the genome files and add them to the channel
    // in order that the files are correctly mounted when using docker
    // is also the opportunity to check that they exists
    input_file.eachLine { line ->
        def genome_file = file(line.split('\t')[1])

        if (!genome_file.exists()) {
            exit 1, "ERROR: Please check input genomes -> Genome file does not exist!\n${genome_file}"
        }

        genome_files.add(file(genome_file))
    }


    // Getting the extension of the first genome file

    extension_patern = ~/.*(\.[a-yA-Y]+)(\.gz)?$/ // A-Y to exclude Z to not cacth gz if exists.
    genome_extension = (genome_files[0] =~ extension_patern)[0][1].toLowerCase()
    annotation_extensions = params.annotation_extensions.split(';')
    fasta_extensions = params.fasta_extensions.split(';')

    if (annotation_extensions.contains(genome_extension)){
        meta.file_type = "annotation"
    } else if (fasta_extensions.contains(genome_extension) ) {
        meta.file_type = "fasta"
    }
    else {
        exit 1, """
        ERROR: Please check input genomes -> Genome file (${genome_files[0]}) does have an unexpected extension: $genome_extension
        Possible value for annotation files: $annotation_extensions\
        Fasta files: $fasta_extensions
        """
    }
    input_meta = [ meta, input_file]

    return input_meta
}


    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    // INPUT_CHECK (
    //     file(params.input)
    // )
    // ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // //
    // // MODULE: Run FastQC
    // //
    // FASTQC (
    //     INPUT_CHECK.out.reads
    // )
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

    // //
    // // MODULE: MultiQC
    // //
    // workflow_summary    = WorkflowPangbank.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // methods_description    = WorkflowPangbank.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    // ch_methods_description = Channel.value(methods_description)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList()
    // )
    // multiqc_report = MULTIQC.out.report.toList()
// }



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    // TODO try too use hook url to send message on mattermost https://developers.mattermost.com/integrate/webhooks/incoming/
    // if (params.hook_url) {
    //     NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    // }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
