process qc_report {

  publishDir "${params.outdir}/qc-report", mode: 'copy'

  input:
    path qc_files_ch                 

  output:
    path 'qc-report.html'          ,  emit: qc_panel_html
    path 'qc-report',  type: 'dir' , optional: true, emit: qc_panel_assets

  script:
    """
    qc_files_path=\$(realpath ${qc_files_ch})
    panel_csv=\$(realpath ${params.panel})
      ln -s ${baseDir}/scripts/qc.qmd qc.qmd
   # panel_csv = ${params.panel} 
    quarto render qc.qmd \
        --execute-dir . \
        --to html \
        --output qc-report.html \
        --output-dir . \
        -P file_list="\${qc_files_path}" \
        -P panel_data="\${panel_csv}"
   """
}