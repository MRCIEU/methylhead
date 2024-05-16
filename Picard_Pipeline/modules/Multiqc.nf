process Multiqc {
        
        input: 
        
        path 'logs/*'
          
        publishDir "${params.outdir}/Multiqc" , mode: 'copy'
      
       
       output:
       
        path "multiqc_report.html", emit: html
        path "multiqc_data", emit: data
       
               
        script:
        """
        mkdir -p ${baseDir}/${params.outdir}/Multiqc
        mkdir -p ${baseDir}/${params.outdir}/Multiqc/multiqc_data
        multiqc .       
        """

}