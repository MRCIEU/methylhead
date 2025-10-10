process multiqc {
  
  input: 
    path 'multiqc_files/*'
          
  output: 
    path "multiqc_report.html"                           , emit: html
    path "multiqc_data"                                  , emit: data
    path "multiqc_data/cutadapt_filtered_reads_plot.txt" , emit: reads
    path "multiqc_data/picard_hsmetrics_table.txt"       , emit: hs
        
  script:
  """      
   multiqc .       
  """
}