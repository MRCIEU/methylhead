#!/usr/bin/env nextflow

/*
 * BISMARK
 */

include { Bismark_pipeline } from './workflows/Bismark_pipeline' 

workflow {
    Bismark_pipeline(params.reads,"${params.outdir}")
}
