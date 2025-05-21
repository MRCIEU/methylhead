#!/usr/bin/env nextflow
      
include { pipeline } from './workflows/pipeline'
    
workflow {
    reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
    outdir = params.outdir
    pipeline(reads, outdir) 
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
