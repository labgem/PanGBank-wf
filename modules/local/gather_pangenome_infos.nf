process GATHER_PANGENOME_INFO {
    label 'process_single'

    // reuse ppanggolin env as it as already been downloaded and used
    conda "bioconda::ppanggolin=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.2.3--hbcbf7aa_0' :
        'biocontainers/ppanggolin:2.2.3--hbcbf7aa_0' }"

    input:
    path "pangenome_infos/*"
    path "genome_stats_summaries/*"


    output:
    path "pangenome_summary.tsv", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gather_pangenome_infos.py --yaml_info_dir pangenome_infos/ --output pangenome_summary.tsv  --yaml_genome_stats_dir genome_stats_summaries/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}

