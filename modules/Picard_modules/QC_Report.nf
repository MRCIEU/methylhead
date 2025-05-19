process QC_Report {

    publishDir "${params.outdir}/QC_Report", mode: 'copy'

    input:
        path qc_files_ch                 

    output:
        path 'qc-report_full.html'  ,               emit: qc_full_html
        path 'qc-report_panel.html' ,               emit: qc_panel_html
        path 'qc-report_full_files' ,               type: 'dir', optional: true, emit: qc_full_assets
        path 'qc-report_panel_files',               type: 'dir', optional: true, emit: qc_panel_assets

    script:
    """
    qc_files_path=\$(realpath ${qc_files_ch})
    panel_csv=\$(realpath ${params.panel_qc})
    
    ln -s ${baseDir}/bin/qc.qmd qc.qmd
    quarto render qc.qmd \
        --execute-dir . \
        --to html \
        --output qc-report_full.html \
        --output-dir . \
        -P file_list="\${qc_files_path}" \
        -P dataset="full"
    
    # panel_csv = ${params.panel_qc} 
    quarto render qc.qmd \
        --execute-dir . \
        --to html \
        --output qc-report_panel.html \
        --output-dir . \
        -P file_list="\${qc_files_path}" \
        -P dataset="panel" \
        -P panel_data="\${panel_csv}"
    """
}