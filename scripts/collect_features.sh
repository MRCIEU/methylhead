#!/bin/bash

set -euo pipefail

SEP=""

while getopts "t" opt; do
  case $opt in
    t)
      SEP="-t"
      ;;
    \?)
      echo "Usage: $0 [-t] arg1 arg2 arg3"
      exit 1
      ;;
  esac
done

shift $((OPTIND - 1))

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 [-t] <files.csv> <columns> <output.csv>"
  exit 1
fi

INPUT="$1"
COLS="$2"
OUTPUT="$3"

FREQ_COLS=$(seq -s, $(awk -F, '{print NF}' <<< "$COLS"))

echo "collect_features.sh---------------------------------------------"

TMPFILE=$(mktemp /tmp/tmpfile.XXXXXX)

cut -d , -f 2 ${INPUT} | tail -n +2 > ${TMPFILE}

echo "$TMPFILE $COLS $FREQ_COLS $OUTPUT"

csvtk concat -X ${TMPFILE} \
    | csvtk cut $SEP -f ${COLS} \
    | csvtk freq $SEP -f ${FREQ_COLS}  \
    > ${OUTPUT}

#rm ${TMPFILE}
