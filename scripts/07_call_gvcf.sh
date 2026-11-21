#!/usr/bin/env bash
# Step 7: per-sample variant calling in GVCF mode, restricted to exome targets (+100 bp padding)
set -euo pipefail
source config.sh
FA=$REF/chr$CHROM.fa
BED=$REF/targets.chr$CHROM.bed
mkdir -p $WORK/gvcf $WORK/tmp

for S in $SAMPLES; do
  if [ -s $WORK/gvcf/$S.g.vcf.gz.tbi ]; then echo "== $S already done, skipping =="; continue; fi
  echo "== HaplotypeCaller on $S (started $(date +%H:%M)) =="
  gatk --java-options "-Xmx6g" HaplotypeCaller \
    -R $FA -I $WORK/bam/$S.dedup.bam \
    -L $BED -ip 100 \
    -ERC GVCF \
    --tmp-dir $WORK/tmp \
    -O $WORK/gvcf/$S.g.vcf.gz 2> $WORK/gvcf/$S.hc.log
  echo "   finished $(date +%H:%M)"
done
ls -lh $WORK/gvcf/*.g.vcf.gz
