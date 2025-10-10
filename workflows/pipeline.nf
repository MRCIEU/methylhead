#!/usr/bin/env nextflow

include { fastqc } from '../modules/fastqc'
include { merge_lanes } from '../modules/merge_lanes'
include { trim_galore } from '../modules/trim_galore'
include { panel_intervals } from '../modules/panel_intervals'
include { bwa_align } from '../modules/bwa_align'
include { sambamba } from '../modules/sambamba' 
include { sort_bam } from '../modules/sort_bam'
include { mark_duplicated } from '../modules/mark_duplicated' 
include { collect_hs_metrics } from '../modules/collect_hs_metrics' 
include { collect_mm_metrics } from '../modules/collect_mm_metrics' 
include { mbias } from '../modules/mbias'
include { extract_cpg_bedgraph } from '../modules/extract_cpg_bedgraph'
include { reformat_bedgraph } from '../modules/reformat_bedgraph'
include { samtools_stats } from '../modules/samtools_stats'
include { extract_cpg_methylation } from '../modules/extract_cpg_methylation'
include { bsmap_align } from '../modules/bsmap_align'
include { camda } from '../modules/camda'
include { methylation_matrix } from '../modules/methylation_matrix'
include { illumina_matrix } from '../modules/illumina_matrix'
include { estimate_cell_counts } from '../modules/estimate_cell_counts'
include { dna_methylation_scores } from '../modules/dna_methylation_scores'
include { camda_matrix } from '../modules/camda_matrix'
include { multiqc } from '../modules/multiqc'
include { qc_report } from '../modules/qc_report'
include { test_associations } from '../modules/test_associations' 
 
workflow pipeline {
    
    take:
    samplesheet
    assembly
    genome_fasta
    panel    
    phenotype    
    models
    cell_reference
    samtools_path    
    outdir
       
    main:

    def required_files = [samplesheet, genome_fasta, panel, phenotype, models, cell_reference]    
    required_files.each { fname ->    
        if (!file(fname).exists()) error "Required file not found: ${fname}"
    }

    read_pairs_ch = Channel
        .fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def sample_id = row.sample_id
            def read1 = file(row.read1, checkIfExists: true)
            def read2 = file(row.read2, checkIfExists: true)
            return tuple(sample_id, read1, read2)
        }
        .groupTuple()
        .map { sample_id, reads_list1, reads_list2 ->
            tuple(sample_id, reads_list1, reads_list2)
        }
	
    merge_lanes(read_pairs_ch)

    clean_files_ch = merge_lanes.out
       .map { sample_id, read1, read2 ->
           [sample_id, [read1, read2]]
        }
	.filter { sample_id, reads -> 
           def fileSize = reads.collect { it.size() / (1024 * 1024) } 
           fileSize[0] >= 0.500 && fileSize[1] >= 0.500
        }
        
    fastqc(clean_files_ch)
    
    trim_galore(clean_files_ch)

    panel_intervals(panel, genome_fasta, "${genome_fasta}.dict") 

    bwa_align(trim_galore.out.fq, genome_fasta)

    sambamba(bwa_align.out.bam, genome_fasta)
	
    sort_bam(sambamba.out.bam)

    mark_duplicated(sort_bam.out.bam)
	
    samtools_stats(
	bwa_align.out.bam,
	mark_duplicated.out.bam,
	sort_bam.out.bam,
	sort_bam.out.bai,
	panel_intervals.out.bed)

    collect_hs_metrics(mark_duplicated.out.bam, genome_fasta, panel_intervals.out.list)

    collect_mm_metrics(mark_duplicated.out.bam, genome_fasta)

    Channel.empty()
	.mix( fastqc.out )             
	.mix( trim_galore.out )        
	.mix( mark_duplicated.out )    
	.mix( collect_hs_metrics.out ) 
	.mix( collect_mm_metrics.out )     
	.mix( samtools_stats.out )                                          
	.map { sample_id, files -> files }
	.collect()
	.set { multiqc_files_ch }
    multiqc(multiqc_files_ch)

    mbias(mark_duplicated.out.bam, genome_fasta)

    extract_cpg_bedgraph(mark_duplicated.out.bam, genome_fasta)
	
    reformat_bedgraph(extract_cpg_bedgraph.out)	

    extract_cpg_methylation(mark_duplicated.out.bam, genome_fasta)

    meth_files = extract_cpg_methylation.out.cpg
	.map { file -> file.toString() }
	.collectFile(name:"meth_files.csv",newLine:true)
    methylation_matrix(meth_files, assembly)

    estimate_cell_counts(methylation_matrix.out.meth, cell_reference)
	
    illumina_matrix(methylation_matrix.out.meth, assembly)
	
    dna_methylation_scores(illumina_matrix.out.meth, assembly)

    bsmap_align(trim_galore.out.fq, genome_fasta)
	
    camda(bsmap_align.out.bam,genome_fasta,samtools_path)
		
    camda_files = camda.out.camda
	.map { file -> file.toString() }
	.collectFile(name:"camda_files.csv",newLine:true)
    camda_matrix(camda_files, assembly)
		
    output_files = multiqc.out.reads
	.concat(estimate_cell_counts.out.counts)
	.concat(methylation_matrix.out.coverage)
	.concat(methylation_matrix.out.meth)
	.concat(illumina_matrix.out.meth)
	.concat(camda_matrix.out.scores)
	.concat(dna_methylation_scores.out.scores)
	.concat(multiqc.out.hs)
	.concat(samtools_stats.out.read_counts_file.map { id, f -> f }) 
	.map { file -> file.toAbsolutePath().toString() }
	.collectFile(name: "output_files.csv", newLine: true)

    qc_report(output_files, panel)
	
    test_associations(output_files, phenotype, models)
}


