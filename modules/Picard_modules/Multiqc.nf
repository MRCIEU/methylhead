process Multiqc {
        
        input: 
        
        path 'multiqc_files/*'
          
        publishDir "${params.outdir}/Multiqc" , mode: 'copy'
      
       
       output:
       
        path "multiqc_report.html", emit: html
        path "multiqc_data", emit: data
       
               
        script:
        """      
        multiqc .       
        """

}