#!/usr/bin/env bash
set -euo pipefail

# Usage: bash reference-create.sh hg19.fa regions.bed output.fa

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <reference_fasta> <regions_bed> <output_fasta>"
    exit 1
fi

REF="$1"            # e.g., hg19.fa
BED="$2"            # e.g., regions.bed
OUTFA="$3"          # e.g., reference.fa

############################################################################
# 1. Make sure coordinates are non-negative and write a temporary BED
############################################################################
awk '
{
  start = $2
  if (start < 0) start = 0
  end   = $3
  print $1 "\t" start "\t" end
}' "$BED" > bed_extended.bed

############################################################################
# 2. Extract those regions from reference FASTA and create a multi-line FASTA
############################################################################
bedtools getfasta -fi "$REF" -bed bed_extended.bed -fo regions_tmp.fa

############################################################################
# 3. Flatten each chromosome’s sequence to a single line (for BWA-meth)
############################################################################
awk '
/^>/ {
  split($0, arr, "[:]")
  chr = substr(arr[1], 2)
  if (!seq[chr]) order[++n] = chr
  curr_chr = chr
  next
}
{
  seq[curr_chr] = seq[curr_chr] $0
}
END {
  for (i = 1; i <= n; i++) {
    print ">" order[i]
    s = seq[order[i]]
    for (j = 1; j <= length(s); j += 60)
      print substr(s, j, 60)
  }
}
' regions_tmp.fa > "$OUTFA"

# Clean up
rm regions_tmp.fa bed_extended.bed

############################################################################
# 4. Create samtools index (.fai) and Picard dictionary (.dict)
############################################################################
samtools faidx "$OUTFA"           > /dev/null
samtools dict "$OUTFA"            > "${OUTFA%.fa}.dict"

############################################################################
# 5. (Optional) quick demo BED: chromosome lengths from .fai, 100bp-10bp windows
############################################################################
awk 'NR>1 { print $1 "\t" ($2-100) "\t" ($2-10) }' "${OUTFA}.fai" > test-target.bed

############################################################################
# 6. Build BWA-Meth index
############################################################################
bwameth.py index "$OUTFA"

echo "Reference built and indexed: $OUTFA"
echo "Index files created: ${OUTFA}.fai, ${OUTFA%.fa}.dict"