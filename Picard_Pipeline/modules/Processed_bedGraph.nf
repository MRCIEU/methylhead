process Processed_bedGraph {
      
    input:
    path bedGraph

    publishDir "${params.outdir}/Processed_bedGraph/", mode: 'copy'

    output:
    path "${bedGraph}_processed.bedgraph", emit: Processed_bedGraph

    script:
    """
    awk 'BEGIN {FS=OFS=\"\t\"} NR == 1 {print \$0} NR > 1 {print \$1,\$2,\$3,((\$5/(\$5+\$6)*100)+0),\$5,\$6;}' OFMT="%.2f" ${bedGraph} > "${bedGraph}_processed.bedgraph"
    """
}
