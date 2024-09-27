process CAMDA {
     
     input:
     tuple val(sample_id), path(bam_camda)
     
     publishDir "${params.outdir}/CAMDA" , mode: 'copy'
     
     output:
     
     tuple val(sample_id), path("${sample_id}_CpG_CAMDA.tsv") , emit: camda
   
     shell:
     """ 
    python ${baseDir}/bin/CAMDA.py CAMDA ${bam_camda} ${params.genome_folder}/hg19.fa -o ${sample_id} -w ${sample_id} -s ${params.samtools_path} -X CG
     """
}
