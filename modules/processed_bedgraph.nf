process processed_bedgraph {
  publishDir "${params.outdir}/processed-bedgraph-files/", mode: 'copy'

  input:
    tuple val(sample_id), path (bedgraph)

  output:
    tuple val(sample_id), path ("${sample_id}_processed.bedgraph") , emit: processed_bedgraph

  script:
    """
    awk 'BEGIN {FS=OFS=\"\t\"} NR == 1 {print \$0} NR > 1 {print \$1,\$2,\$3,((\$5/(\$5+\$6)*100)+0),\$5,\$6;}' OFMT="%.2f" ${bedgraph} > "${sample_id}_processed.bedgraph"
    """
}
