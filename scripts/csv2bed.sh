#!/bin/bash

PANEL=$1
BED=$2

set -euo pipefail

awk 'BEGIN {
    FS = ","
}

NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "chr") { nchr = i }
        if ($i == "start") { nstart = i }
        if ($i == "end") { nend = i }
    }
    next
}

{
    chr = $nchr
    start = $nstart
    end = $nend
    if (chr !~ /^chr/) {
        chr = "chr" chr
    }
    if (start == end) {
        end = end + 1
    }
    printf "%s\t%d\t%d\n", chr, start, end;
}' ${PANEL} |
    sort -k1,1 -k2,2n > ${BED}
