#!/usr/bin/env nextflow
      
include { pipeline } from './workflows/pipeline'
    
workflow {
    samplesheet = params.samplesheet
    assembly = params.assembly
    genome_fasta = params.genome_fasta
    panel = params.panel
    phenotype = params.phenotype
    models = params.models
    cell_reference = params.cell_reference
    samtools_path = params.samtools_path
    outdir = params.outdir

    pipeline(
        samplesheet,
	assembly,
	genome_fasta,
	panel,
	phenotype,
	models,
	cell_reference,
	samtools_path,
	outdir
    ) 
}   
    
workflow.onComplete {

    def msg = """\

        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()
}
