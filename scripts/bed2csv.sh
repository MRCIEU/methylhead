#!/bin/bash

IN=$1
OUT=$2

echo "chr,start,end" > $OUT
cut -f1-3 $IN | tr '\t' ',' >> $OUT
