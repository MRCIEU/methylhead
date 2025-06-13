process interval_file {
  
  publishDir "${params.outdir}/interval-file/", mode: 'copy'

  input:
    path panel
    path params.genome_folder

  output:
    path "interval_file", emit: interval_file

  script:
    """
    set -euo pipefail

    awk -F, 'NR>1 && \$2!=\$3 {printf "%s\\t%d\\t%d\\n",\$1,int(\$2),int(\$3)}' "${panel}" \\
      | sort -k1,1 -k2,2n > panel.bed

    if [[ ! -s panel.bed ]]; then
        echo "panel.bed is empty or missing!"
        exit 1
    fi
    picard BedToIntervalList \
        I=panel.bed \
        O=interval_file \
        SD=${params.genome_folder}
    """
}
