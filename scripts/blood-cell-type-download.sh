#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <output-folder>"
  exit 1
fi

OUTDIR="$1"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

URLS=(
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652277/suppl/GSM5652277%5FBlood%2DT%2DCD3%2DZ000000TV.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652278/suppl/GSM5652278%5FBlood%2DT%2DCD3%2DZ000000UP.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652279/suppl/GSM5652279%5FBlood%2DT%2DCD4%2DZ000000TT.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652280/suppl/GSM5652280%5FBlood%2DT%2DCD4%2DZ000000U7.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652281/suppl/GSM5652281%5FBlood%2DT%2DCD4%2DZ000000UM.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652282/suppl/GSM5652282%5FBlood%2DT%2DCD8%2DZ000000TR.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652283/suppl/GSM5652283%5FBlood%2DT%2DCD8%2DZ000000U5.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652284/suppl/GSM5652284%5FBlood%2DT%2DCD8%2DZ000000UK.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652285/suppl/GSM5652285%5FBlood%2DT%2DCenMem%2DCD4%2DZ00000417.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652286/suppl/GSM5652286%5FBlood%2DT%2DCenMem%2DCD4%2DZ0000041D.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652287/suppl/GSM5652287%5FBlood%2DT%2DCenMem%2DCD4%2DZ0000041N.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652288/suppl/GSM5652288%5FBlood%2DT%2DEff%2DCD8%2DZ00000419.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652289/suppl/GSM5652289%5FBlood%2DT%2DEff%2DCD8%2DZ0000041F.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652290/suppl/GSM5652290%5FBlood%2DT%2DEff%2DCD8%2DZ0000041Q.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652291/suppl/GSM5652291%5FBlood%2DT%2DEffMem%2DCD4%2DZ00000416.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652292/suppl/GSM5652292%5FBlood%2DT%2DEffMem%2DCD4%2DZ0000041C.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652293/suppl/GSM5652293%5FBlood%2DT%2DEffMem%2DCD4%2DZ0000041M.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652294/suppl/GSM5652294%5FBlood%2DT%2DEffMem%2DCD8%2DZ0000041A.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652295/suppl/GSM5652295%5FBlood%2DT%2DEffMem%2DCD8%2DZ0000041G.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652296/suppl/GSM5652296%5FBlood%2DT%2DNaive%2DCD4%2DZ0000041E.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652297/suppl/GSM5652297%5FBlood%2DT%2DNaive%2DCD8%2DZ0000041B.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652298/suppl/GSM5652298%5FBlood%2DT%2DNaive%2DCD8%2DZ0000041H.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652299/suppl/GSM5652299%5FBlood%2DNK%2DZ000000TM.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652300/suppl/GSM5652300%5FBlood%2DNK%2DZ000000U1.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652301/suppl/GSM5652301%5FBlood%2DNK%2DZ000000UF.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652302/suppl/GSM5652302%5FBlood%2DMonocytes%2DZ000000TP.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652303/suppl/GSM5652303%5FBlood%2DMonocytes%2DZ000000U3.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652304/suppl/GSM5652304%5FBlood%2DMonocytes%2DZ000000UH.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652313/suppl/GSM5652313%5FBlood%2DGranulocytes%2DZ000000TZ.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652314/suppl/GSM5652314%5FBlood%2DGranulocytes%2DZ000000UD.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652315/suppl/GSM5652315%5FBlood%2DGranulocytes%2DZ000000UT.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652316/suppl/GSM5652316%5FBlood%2DB%2DZ000000TX.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652317/suppl/GSM5652317%5FBlood%2DB%2DZ000000UB.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652318/suppl/GSM5652318%5FBlood%2DB%2DZ000000UR.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652319/suppl/GSM5652319%5FBlood%2DB%2DMem%2DZ0000041J.beta"
"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/GSM5652320/suppl/GSM5652320%5FBlood%2DB%2DMem%2DZ0000041K.beta"
)

for url in "${URLS[@]}"; do
    wget -c "$url"
done
