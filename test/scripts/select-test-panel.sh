#!/usr/bin/env bash
set -euo pipefail

# USAGE: bash script.sh <BAM_DIR> <BED1> <BED2> <OUT_BED>
# Example: bash script.sh ./bam-files bed1.bed bed2.bed output.bed

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 <BAM_DIR> <BED1> <BED2> <OUT_BED>"
  exit 1
fi

BAM_DIR="$1"
BED1="$2"
BED2="$3"
OUTBED="$4"

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

# 1. Convert each BAM to BED
for bam in "$BAM_DIR"/*.bam; do
  base=$(basename "${bam%.bam}")
  bedtools genomecov -ibam "$bam" -bga > "$TMPDIR/${base}.bed"
done

# 2. Merge BED1 and BED2, remove duplicates
cat "$BED1" "$BED2" | sort -k1,1 -k2,2n -k3,3n | uniq > "$TMPDIR/merged_targets.bed"

# 3. Intersect each BED with merged targets
for bed in "$TMPDIR"/*.bed; do
  base=$(basename "${bed%.bed}")
  bedtools intersect -a "$bed" -b "$TMPDIR/merged_targets.bed" > "$TMPDIR/${base}_covered.bed"
done

# 4. Combine, sort, deduplicate all covered BEDs
cat "$TMPDIR"/*_covered.bed \
  | awk '$2 != $3' \
  | sort -k1,1 -k2,2n -k3,3n \
  | uniq > "$TMPDIR/final-regions.bed"

# 5. Select regions with depth=0 and length>0
awk '$4 == 0 && $2 != $3 {print $1, $2, $3}' "$TMPDIR/final-regions.bed" > "$TMPDIR/final_clean.bed"

# 6. Group by (chr, end), keep minimum start
awk '
{
  key = $1" "$3
  if (!(key in min) || $2 < min[key]) min[key]=$2
}
END {
  for (k in min) {
    split(k, arr, " ")
    print arr[1], min[k], arr[2]
  }
}
' "$TMPDIR/final_clean.bed" | sort -k1,1 -k2,2n > "$TMPDIR/grouped_clean.bed"

# 7. Output as tab-delimited
awk '{$1=$1}1' OFS='\t' "$TMPDIR/grouped_clean.bed" > "$OUTBED"

echo "Output: $OUTBED"
