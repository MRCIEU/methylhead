process Multiqc {
        
        input: 

        path 'multiqc_files/*'
          
        publishDir "${params.outdir}/Multiqc" , mode: 'copy'
          
       output:
       
        path "multiqc_report.html", emit: html
        path "multiqc_data", emit: data
        path "multiqc_data/cutadapt_filtered_reads_plot.txt", emit: reads
           
        script:
        """      
        multiqc .       
        """
}