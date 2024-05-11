#!/usr/bin/env nextflow

/*
 * PICARD
 */

include { Picard_pipeline  } from './workflows/Picard_pipeline' 

workflow {
    Picard_pipeline(params.reads,"${params.outdir}")
}
