process merge_lanes {

    input:
    tuple val(sample_id), path(reads1), path(reads2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz"), emit: pairs

    script:
    """
    cat ${reads1} > ${sample_id}_R1.fastq.gz
    cat ${reads2} > ${sample_id}_R2.fastq.gz
    """
}

