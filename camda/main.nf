#!/usr/bin/env nextflow

/*
 * CAMDA
 */

include { camda_pipeline } from './workflows/camda_pipeline' 

workflow {
    camda_pipeline(params.reads, "${params.outdir}")
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
