#!/usr/bin/env bash
set -euo pipefail

# Check for minimum required arguments: at least one BAM and one BED file
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <bam-files/*.bam> <*.bed>"
    exit 1
fi

bam_files=()
bed_files=()
for arg in "$@"; do
    if [[ "$arg" == *.bam ]]; then
        bam_files+=("$arg")
    elif [[ "$arg" == *.bed ]]; then
        bed_files+=("$arg")
    fi
done

if [[ ${#bam_files[@]} -eq 0 ]]; then
    echo "No BAM files provided!"
    exit 1
fi
if [[ ${#bed_files[@]} -eq 0 ]]; then
    echo "No BED files provided!"
    exit 1
fi

tmp_combined="tmp_combined_$$.bed"

# Create a combined BED from all BAMs without any intermediate folders
for bam in "${bam_files[@]}"; do
    bedtools bamtobed -i "$bam" | cut -f1-3
done | sort -k1,1 -k2,2n | bedtools merge > "$tmp_combined"

# Merge the combined BAM-derived BEDs with additional BEDs
cat "$tmp_combined" "${bed_files[@]}" | sort -k1,1 -k2,2n | bedtools merge > reference-regions.bed

rm -f "$tmp_combined"

echo "Created reference-regions.bed"
