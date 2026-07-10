process MERGE_GTDB_SPLIT_SPECIES {
    tag "$genome_list.baseName"
    label 'process_medium'

    container "ghcr.io/labgem/pangbank-wf:merge-split2"

    input:
    tuple val(meta), path(genome_list)
    path genome_to_fna_paths
    val ani_threshold
    val af_threshold
    val species_ani_stat

    output:
    path "${genome_list.baseName}.clusters", emit: split_clusters
    path "${genome_list.baseName}.genomes.clusters", emit: genome_clusters
    path "${genome_list.baseName}.merge_summary.tsv", emit: merge_summary
    path "${genome_list.baseName}.species_ani.tsv", emit: species_ani
    path "${genome_list.baseName}.species_pair_summary.tsv", emit: species_pair_summary
    path "versions.yml", emit: versions

    script:
    skani_ani_threshold = ani_threshold - 5
    """

    # TODO: the logic could be handle in a python script using pyskani to simplify the process.

    # Get unique species from genome list
    cut -f1 ${genome_list} | sort -u | grep -v '^\$' > species.list

    mkdir -p sketches

    # For each species, extract genome fna paths and create per-species sketches
    while IFS= read -r species; do
        safe_species=\$(echo "\$species" | tr ' /' '__')
        awk -F'\\t' -v sp="\$species" \\
            'NR==FNR { if (\$1==sp) ids[\$2]=1; next } \$1 in ids { print \$2 }' \\
            ${genome_list} ${genome_to_fna_paths} > "\${safe_species}.fna_paths.list"

        skani sketch -t ${task.cpus} \\
            -l "\${safe_species}.fna_paths.list" \\
            -o "sketches/\${safe_species}" \\
            --separate-sketches
        ls "sketches/\${safe_species}"/*.sketch > "\${safe_species}.sketch.list"
    done < species.list

    # Run skani dist for each inter-species pair
    mapfile -t species_array < species.list
    n_species=\${#species_array[@]}
    header_written=0
    for (( i=0; i<n_species; i++ )); do
        for (( j=i+1; j<n_species; j++ )); do
            sp_a=\$(echo "\${species_array[i]}" | tr ' /' '__')
            sp_b=\$(echo "\${species_array[j]}" | tr ' /' '__')
            echo "Calculating distances between \${sp_a} and \${sp_b}"
            if [ "\$header_written" -eq 0 ]; then
                skani dist -s $skani_ani_threshold -t ${task.cpus} --ql "\${sp_a}.sketch.list" --rl "\${sp_b}.sketch.list" \\
                    > ${genome_list.baseName}.dist.tsv
                header_written=1
            else
                skani dist -s $skani_ani_threshold -t ${task.cpus} --ql "\${sp_a}.sketch.list" --rl "\${sp_b}.sketch.list" \\
                    | tail -n +2 >> ${genome_list.baseName}.dist.tsv
            fi
        done
    done

    merge_gtdb_splits.py --genome-list ${genome_list} \\
                         --skani-dist ${genome_list.baseName}.dist.tsv \\
                         --genome-fna-paths ${genome_to_fna_paths} \\
                         --ani_threshold ${ani_threshold} \\
                         --af_threshold ${af_threshold} \\
                         --species-ani-stat ${species_ani_stat} \
                         --prefix ${genome_list.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        skani: \$(skani --version 2>&1 | sed 's/^.*skani //; s/ .*\$//')
    END_VERSIONS
    """
}
