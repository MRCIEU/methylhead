#!/usr/bin/env bash
# =============================================================================
#  Prepare hg19 reference genome and optionally send an email on completion.
#
#  Usage examples:
#    With e-mail ─>  bash code.sh -N you@example.com /path/to/ref
#    No e-mail  ──>  bash code.sh /path/to/ref
#
#  Options:
#    -N <email>   Address to notify. If -N is omitted or empty, no mail is sent.
#
#  Requirements:
#    • Bash 4+, wget, apptainer (Singularity), samtools, picard, bwameth.py
#    • A working “mail”/“mailx” command (e.g. mailutils, msmtp, postfix…)
# =============================================================================
set -euo pipefail

# -----------------------------------------------------------------------------
#  ARGUMENT PARSING
# -----------------------------------------------------------------------------
EMAIL=""
while getopts ":N:" opt; do
    case $opt in
        N) EMAIL="$OPTARG" ;;
        \?|:)
            printf "Usage: %s [-N email] [REF_DIR]\n" "$0" >&2
            exit 2
            ;;
    esac
done
shift $((OPTIND-1))                 # positional parameter(s) that remain

# -----------------------------------------------------------------------------
#  SETTINGS
# -----------------------------------------------------------------------------
REF_DIR=${1:-reference}
FA=hg19.fa
URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/${FA}.gz"

SIF_NAME="methylhead-pipeline_wgbs_analysis.sif"
SIF_ORAS="oras://docker.io/onuroztornaci/methylhead-pipeline:wgbs_analysis"

LOGFILE="$(pwd)/$(basename "$0").log" 
SECONDS=0                            # built-in timer for elapsed seconds

# -----------------------------------------------------------------------------
#  PIPELINE
# -----------------------------------------------------------------------------
exec > >(tee -a "$LOGFILE") 2>&1     # log stdout & stderr to file and screen

mkdir -p "$REF_DIR" && cd "$REF_DIR"

rm -f "$SIF_NAME"
apptainer pull "$SIF_NAME" "$SIF_ORAS"

# 1) Download and decompress hg19 if not present
[ -s "$FA" ] || { wget -qO- "$URL" | gunzip -c > "$FA"; }

# 2) FASTA index (.fai)
[ -s "$FA.fai" ] || apptainer exec "$SIF_NAME" samtools faidx "$FA"

# 3) Sequence dictionary (.dict)
DICT="${FA}.dict"
[ -s "$DICT" ] || apptainer exec "$SIF_NAME" picard CreateSequenceDictionary \
        REFERENCE="$FA" OUTPUT="$DICT"

# 4) BWA-meth index
if [ ! -f "${FA}.c2t.sa" ] || [ "$FA" -nt "${FA}.c2t.sa" ]; then
    apptainer exec "$SIF_NAME" bwameth.py index "$FA"
fi

echo "Reference genome prepared in $(pwd)"

# -----------------------------------------------------------------------------
#  CLEANUP & OPTIONAL E-MAIL
# -----------------------------------------------------------------------------
cleanup() {
    STATUS=$?
    ELAPSED=$SECONDS
    HUMAN_FMT="$((ELAPSED/3600))h $((ELAPSED%3600/60))m $((ELAPSED%60))s"

    # no email requested
    [[ -z "$EMAIL" ]] && return

    if [ "$STATUS" -eq 0 ]; then
        SUBJECT="[SUCCESS] hg19 prep finished (${HUMAN_FMT})"
        BODY="The hg19 reference genome has been prepared successfully.
Output directory: $(pwd)
Total runtime    : ${HUMAN_FMT}

Regards,
$(hostname)"
    else
        SUBJECT="[FAIL] hg19 prep exited with code ${STATUS} after ${HUMAN_FMT}"
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
