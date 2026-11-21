#!/usr/bin/env bash
# Step 1b: GIAB v4.2.1 truth sets (GRCh37) for all three trio members, chromosome 20 only
set -euo pipefail
source config.sh
mkdir -p $DATA/truth

REL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio
declare -A DIR=( [HG002]=HG002_NA24385_son [HG003]=HG003_NA24149_father [HG004]=HG004_NA24143_mother )

for S in $SAMPLES; do
  echo "== $S truth set =="
  BASE=$REL/${DIR[$S]}/NISTv4.2.1/GRCh37/${S}_GRCh37_1_22_v4.2.1_benchmark
  # variants GIAB is confident about (region query: only chr20 is downloaded)
  bcftools view -r $CHROM -Oz -o $DATA/truth/$S.chr$CHROM.vcf.gz $BASE.vcf.gz
  tabix -p vcf $DATA/truth/$S.chr$CHROM.vcf.gz
  # regions where GIAB is confident; outside these it makes no claim
  curl -fL $BASE"_noinconsistent.bed" | awk -v c=$CHROM '$1==c' > $DATA/truth/$S.chr$CHROM.bed
  echo "$S: $(bcftools view -H $DATA/truth/$S.chr$CHROM.vcf.gz | wc -l) truth variants, $(wc -l < $DATA/truth/$S.chr$CHROM.bed) confident regions"
done
