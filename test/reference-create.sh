#!/usr/bin/env bash
set -euo pipefail

#######################################################################
#  build_reference_from_regions.sh
#
#  What does it do?
#  1. Takes a BED file that represents the *combined* genomic loci from:
#       • the Illumina-450 K methylation-matrix, and
#       • the blood_cell_type.bed panel,
#     then extracts those sequences from the hg19 reference.
#
#  2. Collapses multi-line FASTA entries so each chromosome/contig is on
#     a single line (many aligners and some QC tools prefer this).
#
#  3. Builds samtools and BWA-Meth indices for the trimmed reference.
#
#  Inputs
#  ──────
#    BED   = final_regions.bed   (union described above)
#    REF   = hg19.fa             (full hg19 genome FASTA)
#
#  Outputs
#  ───────
#    reference.fa        – FASTA of only the requested regions
#    reference.fa.fai    – samtools fasta index
#    reference.fa.dict   – Picard sequence dictionary
#    example.bed         – tiny ±-window demo BED built from .fai
#
#  Requirements
#  • awk, bedtools, samtools
#  • bwameth.py  (https://github.com/brentp/bwa-meth)
#######################################################################

BED="final_regions.bed"   # union BED (Illumina-450 K + blood cell-type)
REF="hg19.fa"             # full reference genome
OUTFA="reference.fa"      # trimmed reference to be created

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
# 2. Extract those regions from hg19 and create a multi-line FASTA
############################################################################
bedtools getfasta -fi "$REF" -bed bed_extended.bed -fo regions_tmp.fa

############################################################################
# 3. Flatten each chromosome’s sequence to a single line
#    (makes downstream indexing slightly faster / simpler)
############################################################################
awk '
# Header line: >chr:start-end
/^>/ {
  split($0, arr, "[:]")
  chr = substr(arr[1], 2)         # remove leading ">"
  if (!seq[chr]) order[++n] = chr # remember original order
  curr_chr = chr
  next
}
# Sequence lines
{
  seq[curr_chr] = seq[curr_chr] $0
}
END {
  # Print in original chromosome order, 60-bp wrapped
  for (i = 1; i <= n; i++) {
    print ">" order[i]
    s = seq[order[i]]
    for (j = 1; j <= length(s); j += 60)
      print substr(s, j, 60)
  }
}
' regions_tmp.fa > "$OUTFA"

# House-keeping
rm regions_tmp.fa
rm bed_extended.bed

############################################################################
# 4. Create samtools index (.fai) and Picard dictionary (.dict)
############################################################################
samtools faidx   "$OUTFA"            > /dev/null
samtools dict    "$OUTFA"            > reference.fa.dict

############################################################################
# 5. (Optional) quick demo BED: take chromosome lengths from .fai
#    and build tiny –100 bp/–10 bp windows.  Replace or remove as needed.
############################################################################
awk 'NR>1 { print $1 "\t" ($2-100) "\t" ($2-10) }' reference.fa.fai \
  > example.bed

############################################################################
# 6. Build BWA-Meth index
############################################################################
bwameth.py index "$OUTFA"

echo "Reference built and indexed: $OUTFA"
