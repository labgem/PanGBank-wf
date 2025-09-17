process INDEX_BANK {
    conda '/home/tlemane/work/dev/orgs/LABGeM/PanGBank-wf/envMETAPANG'
    
    input:
        path(index_bank_input_tsv)
        val kmer_size
        val scaled


    output:
        path "bank_index/genome_index.sbt.zip", emit: genome_index_file
        path "bank_index/pangenome_index.sbt.zip", emit: pangenome_index_file
        path "bank_index/index_info.json", emit: index_info_file
        path "stats.txt", emit: stats
        path "versions.yml", emit: versions
        
    script:
    """
    metapang index bank --input ${index_bank_input_tsv} --output bank_index --threads ${task.cpus} --kmer-size ${kmer_size} --scaled ${scaled} 

    SIZEG=\$(stat -c%s "bank_index/genome_index.sbt.zip")
    SIZEP=\$(stat -c%s "bank_index/pangenome_index.sbt.zip")
    SIZEI=\$(stat -c%s "bank_index/index_info.json")

    echo -e "File\tSize" > stats.txt
    echo -e "genome_index.sbt.zip\t\$SIZEG" >> stats.txt
    echo -e "pangenome_index.sbt.zip\t\$SIZEP" >> stats.txt
    echo -e "index_info.json\t\$SIZEI" >> stats.txt

    awk 'BEGIN {OFS="\t"} {if (NR > 1) {\$2 = \$2 / (1024*1024)}; print}' stats.txt > stats_tmp.txt
    mv stats_tmp.txt stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metapang: \$(metapang version --no-name)
    END_VERSIONS
    """
}