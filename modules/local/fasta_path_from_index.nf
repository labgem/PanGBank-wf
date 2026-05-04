// For each genome in an annotation input list, look up its fasta path from a genome_fasta
// index file (genome_name -> fasta_path). Produces the same two output files as ANY2FASTA
// (genome_name_to_fasta.tsv and fasta_to_orginal_path.tsv) without running any conversion.
process FASTA_PATH_FROM_INDEX {
    tag "$meta.species"
    label 'process_single'

    conda "bioconda::ppanggolin=2.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin%3A2.3.0--py312h247cb63_0' :
        'biocontainers/ppanggolin:2.3.0--py312h247cb63_0' }"

    input:
    tuple val(meta), path(annotation_genome_list)
    path genome_fasta_index

    output:
    tuple val(meta), path("genome_name_to_fasta.tsv")  , emit: genome_path_fasta
    tuple val(meta), path("fasta_to_orginal_path.tsv") , emit: fasta_to_orginal_path
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Decompress if .gz, otherwise pass through as-is.
    _dc() { case "\$1" in *.gz) zcat "\$1" ;; *) cat "\$1" ;; esac; }

    awk -F'\\t' '
        NR==FNR { fasta_path[\$1] = \$2; next }
        {
            genome_name = \$1
            ann_path    = \$2
            if (!(genome_name in fasta_path)) {
                missing[genome_name] = 1
            } else {
                print genome_name "\\t" fasta_path[genome_name] >> "genome_name_to_fasta.tsv"
                print fasta_path[genome_name] "\\t" ann_path     >> "fasta_to_orginal_path.tsv"
            }
        }
        END {
            if (length(missing) > 0) {
                msg = "ERROR: " length(missing) " genome(s) not found in --genome_fasta:"
                for (g in missing) msg = msg " " g
                print msg > "/dev/stderr"
                exit 1
            }
        }
    ' <(_dc "${genome_fasta_index}") <(_dc "${annotation_genome_list}")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1)
    END_VERSIONS
    """
}
