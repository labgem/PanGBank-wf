process INDEX_PANGENOME {
    conda "/home/tlemane/work/dev/orgs/LABGeM/PanGBank-wf/envMETAPANG"

    input:
        tuple val(meta), path(pangenome), path(genomes_fof)
        val kmer_size
        val annotation_type

    output:
        path "pangenome_dbg/${meta.species}/*" 
        path "stats.txt", emit: stats
        path "versions.yml", emit: versions

    script:
    def decompressed_genomes = genomes_fof.name.endsWith('.gz') ? "${genomes_fof.baseName}" : genomes_fof
    def genomes = genomes_fof.name != 'NO_FILE' ? "--genomes genomes.txt" : ''
    """
    if [[ "${genomes_fof.name}" == *.gz ]]; then
        gunzip -c ${genomes_fof} > ${decompressed_genomes}
    fi
    cut -f2 ${decompressed_genomes} > genomes.txt
    metapang -v trace index pangenome --pangenome ${pangenome} \
                             --name ${meta.species} \
                             --output pangenome_dbg \
                             --kmer-size ${kmer_size} \
                             --annotation-type ${annotation_type} \
                             --threads ${task.cpus} \
                             --tmp ${meta.species}_tmp \
                            ${genomes}

    DBG_FILE=\$(find "pangenome_dbg/${meta.species}" -maxdepth 1 -type f -name "*.dbg" | head -n1)
    ANNODBG_FILE=\$(find "pangenome_dbg/${meta.species}" -maxdepth 1 -type f -name "*.annodbg" | head -n1)

    SIZE_DBG=\$(stat -c%s "\$DBG_FILE")
    SIZE_ANNODBG=\$(stat -c%s "\$ANNODBG_FILE")
    echo -e "Species\tDBG size\tAnnotation size" > stats.txt
    echo -e "${meta.species}\t\$SIZE_DBG\t\$SIZE_ANNODBG" >> stats.txt

    rm -rf ${meta.species}_tmp

    awk 'BEGIN {OFS="\t"} {if (NR > 1) {\$2 = \$2 / (1024*1024); \$3 = \$3 / (1024*1024)}; print}' stats.txt > stats_tmp.txt
    mv stats_tmp.txt stats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metapang: \$(metapang version --no-name)
    END_VERSIONS
    """
}