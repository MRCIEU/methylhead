#!/usr/bin/env bash
set -euo pipefail

# All processing is performed based on the regions defined in the provided BED file (e.g., blood_cell_types_extended.bed).

############################################
#  ARGUMENTS & DEFAULTS
############################################
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <BED_FILE> <OUTPUT_DIR>"
  exit 1
fi

BED_ORIG=$(readlink -f "$1")        # user-supplied BED
OUT_DIR=$(readlink -f "$2")         # where all results will live
FASTQ_DIR=$(pwd)                    # the directory containing *.fastq.gz pairs
THREADS=8                           # CPU cores for bwa / samtools
EXPAND=1000                          # ±bp to expand BED intervals
REF_DIR="${OUT_DIR}/reference"
REF_FASTA="${REF_DIR}/hg19.fa"

mkdir -p "$OUT_DIR" "$REF_DIR"

############################################
#  DOWNLOAD & INDEX hg19 IF NEEDED
############################################
if [[ ! -f "${REF_FASTA}" ]]; then
  echo "[+] Downloading hg19 reference..."
  wget -q -O "${REF_FASTA}.gz" \
       "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
  gunzip "${REF_FASTA}.gz"

  echo "[+] Building BWA index..."
  bwa index -a bwtsw "${REF_FASTA}"

  echo "[+] Building samtools faidx..."
  samtools faidx "${REF_FASTA}"
fi

############################################
#  EXPAND THE BED FILE
############################################
EXP_BED="${OUT_DIR}/regions_expanded.bed"
awk -v OFS="\t" -v EXP="$EXPAND" '
  { s=$2-EXP; if(s<0)s=0; e=$3+EXP; $2=s; $3=e; print }
' "$BED_ORIG" > "$EXP_BED"
echo "[+] BED expanded -> $EXP_BED"

############################################
#  LOOP OVER FASTQ PAIRS
############################################
for fq1 in "${FASTQ_DIR}"/*_R1*.fastq.gz; do
  [[ -e "$fq1" ]] || { echo "(!) No FASTQ pairs found"; exit 1; }

  fq2="${fq1/_R1/_R2}"
  sample=$(basename "$fq1" | sed 's/_R1.*\.fastq\.gz//')

  if [[ ! -f "$fq2" ]]; then
    echo "(!) Missing R2 for $sample — skipping"
    continue
  fi

  echo "========================================"
  echo "[*] Aligning and filtering: $sample"

  # 1) ALIGN → coordinate-sorted BAM
  bam_sorted="${OUT_DIR}/${sample}.sorted.bam"
  bwa mem -t "$THREADS" "$REF_FASTA" "$fq1" "$fq2" |
    samtools sort -@ "$THREADS" -o "$bam_sorted" -
  samtools index "$bam_sorted"

  # 2) FILTER BY BED REGIONS
  bam_filt="${OUT_DIR}/${sample}.filtered.bam"
  samtools view -@ "$THREADS" -b -L "$EXP_BED" -o "$bam_filt" "$bam_sorted"

  # 3) BAM → FASTQ (name-sorted first)
  fq_out1="${OUT_DIR}/${sample}_1.fastq.gz"
  fq_out2="${OUT_DIR}/${sample}_2.fastq.gz"
  samtools sort -@ "$THREADS" -n -o - "$bam_filt" |
    samtools fastq -@ "$THREADS" -1 "$fq_out1" -2 "$fq_out2" -

  # 4) SYNCHRONISE PAIRS (pure Bash + awk)
  zcat "$fq_out1" | awk 'NR%4==1{gsub(/^@/,"");print}' > "${OUT_DIR}/${sample}.ids1"
  zcat "$fq_out2" | awk 'NR%4==1{gsub(/^@/,"");print}' > "${OUT_DIR}/${sample}.ids2"
  sort "${OUT_DIR}/${sample}.ids1" -o "${OUT_DIR}/${sample}.ids1"
  sort "${OUT_DIR}/${sample}.ids2" -o "${OUT_DIR}/${sample}.ids2"
  comm -12 "${OUT_DIR}/${sample}.ids1" "${OUT_DIR}/${sample}.ids2" \
       > "${OUT_DIR}/${sample}.common"

  for end in 1 2; do
    fq_tmp="${OUT_DIR}/${sample}_${end}.synced.gz"
    zcat "${OUT_DIR}/${sample}_${end}.fastq.gz" | \
      awk -v ids="${OUT_DIR}/${sample}.common" '
        BEGIN{ while((getline<ids)>0) keep[$1]=1 }
        NR%4==1{hdr=$0; id=substr(hdr,2); flag=keep[id]}
        flag{print hdr; getline; print; getline; print; getline; print}
      ' | gzip > "$fq_tmp"
    mv "$fq_tmp" "${OUT_DIR}/${sample}_${end}.fastq.gz"
  done

  rm "${OUT_DIR}/${sample}.ids1" "${OUT_DIR}/${sample}.ids2" \
     "${OUT_DIR}/${sample}.common" "$bam_filt"
  echo "Finished: $sample"
done

echo "========================================"
echo "All samples processed!"
