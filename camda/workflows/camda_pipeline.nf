#!/usr/bin/env nextflow

include { Trim_galore } from '../modules/Trim_galore'
include { BSmap_Aligment } from '../modules/BSmap_Aligment'
include { CAMDA } from '../modules/CAMDA'
include { Camda_matrix } from '../modules/Camda_matrix'

workflow camda_pipeline {
    take:
    reads   
    outdir

    main:
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)    
  Trim_galore(read_pairs_ch)
    trim_ch = Trim_galore.out.fq     
  BSmap_Aligment(trim_ch)
    bam_files_ch = BSmap_Aligment.out.bam
  CAMDA(bam_files_ch)             
     files_ch= CAMDA.out.camda
          .map { file -> file.toString() }
          .collectFile(name:"files.csv",newLine:true)
  Camda_matrix(files_ch)     
}
