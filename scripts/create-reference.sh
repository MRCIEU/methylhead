#!/usr/bin/env bash
# =============================================================================
#  Prepare hg19, hg38, or mm10 reference genome with optional email notification.
#
#  Usage:
#    bash create-reference.sh <genome> <output_dir> [-N <email>]
#
#  Examples:
#    bash create-reference.sh hg19 /your/path/hg19_ref
#    bash create-reference.sh mm10 /anywhere/mm10_ref -N you@email.com
#
#  Arguments:
#    <genome>      One of: hg19 | hg38 | mm10
#    <output_dir>  Output directory (full path, will be created if not exists)
#    -N <email>    (Optional) Email address for notification
# =============================================================================
set -euo pipefail

# Check arguments
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <genome> <output_dir> [-N <email>]"
    echo "       genome: hg19 | hg38 | mm10"
    exit 1
fi

GENOME=$(readlink -f "$1")
OUTDIR=$(readlink -f "$2")
EMAIL=""

if ! [ -f ${GENOME} ]; then
    URL_TEMPLATE="https://hgdownload.soe.ucsc.edu/goldenPath/ASSEMBLY/bigZips/ASSEMBLY.fa.gz"
    URL=${URL_TEMPLATE/ASSEMBLY/$GENOME}
    FA="${GENOME}.fa"
else
    FA=$GENOME
fi

# Optional: email notification
if [[ $# -ge 4 && "$3" == "-N" ]]; then
    EMAIL="$4"
fi

##################################################
## container config
##################################################
SIF_NAME="methylhead-pipeline_wgbs_analysis.sif"
SIF_ORAS="oras://docker.io/onuroztornaci/methylhead-pipeline:wgbs_analysis"

# Download container image and reference genome
rm -f "$SIF_NAME"
apptainer pull "$SIF_NAME" "$SIF_ORAS"


##################################################
## logging config
##################################################
LOGFILE="$(pwd)/$(basename "$0").log"
SECONDS=0

# Log stdout and stderr to file and screen
exec > >(tee -a "$LOGFILE") 2>&1

##################################################
## download and index genome
#################################################

# Create output directory and enter it
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Download and decompress reference genome if not present
[ -s "$FA" ] || { wget -qO- "$URL" | gunzip -c > "$FA"; }

# Create FASTA index (.fai)
[ -s "$FA.fai" ] || apptainer exec "$SIF_NAME" samtools faidx "$FA"

# Create sequence dictionary (.dict)
DICT="${FA}.dict"
[ -s "$DICT" ] || apptainer exec "$SIF_NAME" picard CreateSequenceDictionary \
    REFERENCE="$FA" OUTPUT="$DICT"

# Create BWA-meth index
if [ ! -f "${FA}.c2t.sa" ] || [ "$FA" -nt "${FA}.c2t.sa" ]; then
    apptainer exec "$SIF_NAME" bwameth.py index "$FA"
fi

echo "Reference genome prepared in $(pwd)"

#################################################
# cleanup
#################################################
cleanup() {
    STATUS=$?
    ELAPSED=$SECONDS
    HUMAN_FMT="$((ELAPSED/3600))h $((ELAPSED%3600/60))m $((ELAPSED%60))s"
    rm -f "$SIF_NAME"
    [[ -z "$EMAIL" ]] && return

    if [ "$STATUS" -eq 0 ]; then
        SUBJECT="[SUCCESS] ${GENOME} prep finished (${HUMAN_FMT})"
        BODY="The ${GENOME} reference genome has been prepared successfully.
Output directory: $(pwd)
Total runtime    : ${HUMAN_FMT}

Regards,
$(hostname)"
    else
        SUBJECT="[FAIL] ${GENOME} prep exited with code ${STATUS} after ${HUMAN_FMT}"
        BODY="The pipeline terminated with an error (exit code: ${STATUS}).
Please check the attached log excerpt for details.

Regards,
$(hostname)"
    fi

    {
        printf "%s\n\n" "$BODY"
        echo "----- LOG (last 200 lines) -----"
        tail -n 200 "$LOGFILE"
    } | mail -s "$SUBJECT" "$EMAIL"
}
trap cleanup EXIT
