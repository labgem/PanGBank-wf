process MERGE_GTDB_SPLIT_SPECIES {
    tag "$genome_list.baseName"
    label 'process_medium'

    container "ghcr.io/labgem/pangbank-wf:merge-split2"

    input:
    tuple val(meta), path(genome_list)
    path genome_paths
    val threshold

    output:
    path "${genome_list.baseName}.clusters", emit: split_clusters
    path "${genome_list.baseName}.genomes.clusters", emit: genome_clusters
    path "versions.yml", emit: versions

    script:
    """
    awk 'NR==FNR { ids[\$NF]; next } \$1 in ids { print \$2 }' \
           $genome_list $genome_paths > skani_input.list

    skani triangle -t $task.cpus -l ./skani_input.list --medium -o ${genome_list.baseName}.tsv

    awk '
       NR==FNR {
           path_to_id[\$2] = \$1
           next
       }
       NR==1 {
           for (i=2; i<=NF; i++) {
               if (\$i in path_to_id) {
                   \$i = path_to_id[\$i]
               }
           }
       }
       {
           if (\$1 in path_to_id) {
               \$1 = path_to_id[\$1]
           }
           print
       }
    ' OFS='\\t' $genome_paths ${genome_list.baseName}.tsv > ${genome_list.baseName}.clean.tsv


    merge_gtdb_splits.py --genome-list $genome_list \
                         --skani-triangle ${genome_list.baseName}.clean.tsv \
                         --threshold $threshold \
                         --prefix $genome_list.baseName

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        skani: \$(skani --version 2>&1 | sed 's/^.*skani //; s/ .*\$//')
    END_VERSIONS
    """
}
