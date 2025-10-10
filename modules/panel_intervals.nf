process panel_intervals {
  
  input:
    path panel
    path genome_fasta
    path genome_dict

  output:
    path "interval_file" , emit: list
    path "panel.bed"     , emit: bed 

  script:
    """
    set -euo pipefail

    awk -F, 'NR>1 {
        if (\$2 == \$3)
            printf "%s\\t%d\\t%d\\n", \$1, \$2, \$3+1;
        else
            printf "%s\\t%d\\t%d\\n", \$1, \$2, \$3;
                  }' "${panel}" |
    sort -k1,1 -k2,2n > panel.bed

    if [[ ! -s panel.bed ]]; then
        echo "panel.bed is empty or missing!"
        exit 1
    fi
    picard BedToIntervalList \
        I=panel.bed \
        O=interval_file \
        SD="${genome_fasta}"
    """
}
