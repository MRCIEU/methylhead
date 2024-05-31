process CAMDA {

     
     input:
     tuple val(sample_id), path(bam)
     
     publishDir "${params.outdir}/CAMDA/" , mode: 'copy'
     
     output:
     
     tuple val(sample_id), path("${sample_id}_CpG_MethRatio.tsv")    , emit: MethRatio
     tuple val(sample_id), path("${sample_id}_CpG_CAMDA.tsv")        , emit: CAMDA
     tuple val(sample_id), path("${sample_id}_CpG_MethRatio.wig")    , emit: MethRatio_wig
     tuple val(sample_id), path("${sample_id}_CpG_CAMDA.wig")        , emit: CAMDA_wig

     
     shell:
     """ 
    mkdir -p ${baseDir}/${params.outdir}/CAMDA/     
    python ${workflow.projectDir}/scripts/CAMDA.py CAMDA ${bam} ${params.genome_folder} -o ${sample_id} -w ${sample_id} -s ${params.samtools_path}
     """
}