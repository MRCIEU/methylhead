process download_singularity_image {

     input : 
   
     file 'dnam_cancer_pipeline_latest.sif'

     tag 'Download Singularity Image'

    output:
    file 'dnam_cancer_pipeline_latest.sif'

    script:
    """
    if [ ! -f ${params.sif_file} ]; then
        singularity pull docker://onuroztornaci/dnam_cancer_pipeline:latest
    fi
    """
}

