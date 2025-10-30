#!/bin/bash

PANEL=$1
BED=$2

set -euo pipefail

awk -F, 'NR>1 {
        if ($2 == $3)
            printf "%s\t%d\t%d\n", $1, $2, $3+1;
        else
            printf "%s\t%d\t%d\n", $1, $2, $3;
    }' ${PANEL} |
    sort -k1,1 -k2,2n > ${BED}

