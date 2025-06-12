#!/usr/bin/env bash

######################################################################
#  Prepare reference genome
#
#  Usage:
#    bash scripts/create-reference.sh <genome> <output_dir> 
#
#  Examples:
#    bash create-reference.sh hg19 /your/path/genome/index
#    bash create-reference.sh data/genome.fa /your/path/genome/index
#
#  Arguments:
#    <genome>      The name of a UCSC genome assembly or fasta file
#    <output_dir>  Genome index directory (created if it doesn't exist)
########################################################################

set -euo pipefail

# Check arguments
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <genome> <output_dir>"
    exit 1
fi

GENOME=$(readlink -f "$1")
OUTDIR=$(readlink -f "$2")

## GENOME is not a file that exists, then assume it is a genome assembly name
if ! [ -f ${GENOME} ]; then
    URL_TEMPLATE="https://hgdownload.soe.ucsc.edu/goldenPath/ASSEMBLY/bigZips/ASSEMBLY.fa.gz"
    URL=${URL_TEMPLATE/ASSEMBLY/$GENOME}
    FA="${GENOME}.fa"
else
    FA=$GENOME
fi

# Create output directory and enter it
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Download and decompress reference genome if not present
[ -s "$FA" ] || { wget -qO- "$URL" | gunzip -c > "$FA"; }

# Create FASTA index (.fai)
[ -s "$FA.fai" ] || samtools faidx "$FA"

# Create sequence dictionary (.dict)
DICT="${FA}.dict"
[ -s "$DICT" ] || picard CreateSequenceDictionary REFERENCE="$FA" OUTPUT="$DICT"

# Create BWA-meth index
if [ ! -f "${FA}.c2t.sa" ] || [ "$FA" -nt "${FA}.c2t.sa" ]; then
    bwameth.py index "$FA"
fi

echo "Reference genome prepared in $(pwd)"

