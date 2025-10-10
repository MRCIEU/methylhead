process reformat_bedgraph {

  input:
    tuple val(sample_id), path (bedgraph)

  output:
    tuple val(sample_id), path ("${sample_id}_processed.bedgraph.gz") , emit: processed_bedgraph

  script:
    """
    gunzip -c ${bedgraph} | awk 'BEGIN {FS=OFS=\"\t\"} NR == 1 {print \$0} NR > 1 {print \$1,\$2,\$3,((\$5/(\$5+\$6)*100)+0),\$5,\$6;}' OFMT="%.2f" | gzip > "${sample_id}_processed.bedgraph.gz"
    """
}
