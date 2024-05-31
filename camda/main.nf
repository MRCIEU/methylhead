#!/usr/bin/env nextflow

/*
 * CAMDA
 */

include { camda_pipeline } from './workflows/camda_pipeline' 

workflow {
    camda_pipeline(params.reads,"${params.outdir}")
}