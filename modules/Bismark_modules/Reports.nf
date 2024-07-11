 process Reports {

        input:
        
      	path report_files , from: 'Report_files/*'
                
        publishDir "${params.outdir}/Reports" , mode: 'copy'
      
        output:
        
        path("*.html"), emit: reports
        
        script:
        """
       bismark2report ${report_files}       
        """
}


