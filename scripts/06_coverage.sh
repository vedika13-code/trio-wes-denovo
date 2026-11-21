#!/usr/bin/env bash
# Step 6: depth of coverage over the exome targets (duplicates and low-quality reads excluded)
set -euo pipefail
source config.sh
mkdir -p $RESULTS/align
BED=$REF/targets.chr$CHROM.bed

echo -e "sample\ttarget_bases\tmean_depth\tpct_ge10x\tpct_ge20x\tpct_zero" > $RESULTS/align/06_coverage.tsv
for S in $SAMPLES; do
  echo "== coverage for $S =="
  # -a: report every target base, even those with zero reads
  # -Q 20: count only reads with mapping quality >= 20; -q 20: only bases with quality >= 20
  samtools depth -a -b $BED -Q 20 -q 20 $WORK/bam/$S.dedup.bam \
    | awk -v s=$S '{n++; sum+=$3; if($3>=10)a++; if($3>=20)b++; if($3==0)z++}
        END{printf "%s\t%d\t%.1f\t%.1f%%\t%.1f%%\t%.1f%%\n", s, n, sum/n, 100*a/n, 100*b/n, 100*z/n}' \
    >> $RESULTS/align/06_coverage.tsv
done
cat $RESULTS/align/06_coverage.tsv
