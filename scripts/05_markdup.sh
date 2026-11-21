#!/usr/bin/env bash
# Step 5: flag PCR/optical duplicates so they don't count as independent evidence.
# Uses samtools markdup (GATK MarkDuplicates hit intermittent write failures in the OrbStack VM).
set -euo pipefail
source config.sh
mkdir -p $RESULTS/align $WORK/tmp

echo -e "sample\treads_examined\tduplicate_reads\tpercent_duplication" > $RESULTS/align/05_duplication.tsv
for S in $SAMPLES; do
  echo "== marking duplicates in $S =="
  rm -f $WORK/bam/$S.dedup.bam $WORK/bam/$S.dedup.bai $WORK/bam/$S.dedup.bam.bai
  samtools collate -u -O -T $WORK/tmp/${S}_col $WORK/bam/$S.sorted.bam \
    | samtools fixmate -m -u - - \
    | samtools sort -u -T $WORK/tmp/${S}_srt - \
    | samtools markdup -f $WORK/bam/$S.markdup_stats.txt -T $WORK/tmp/${S}_md - $WORK/bam/$S.dedup.bam
  samtools index $WORK/bam/$S.dedup.bam
  cp $WORK/bam/$S.markdup_stats.txt $RESULTS/align/
  ex=$(awk -F': ' '/^EXAMINED/{print $2}' $WORK/bam/$S.markdup_stats.txt)
  du=$(awk -F': ' '/^DUPLICATE PRIMARY TOTAL/{print $2}' $WORK/bam/$S.markdup_stats.txt)
  awk -v s=$S -v e=$ex -v d=$du 'BEGIN{printf "%s\t%s\t%s\t%.1f%%\n", s, e, d, 100*d/e}' >> $RESULTS/align/05_duplication.tsv
done
cat $RESULTS/align/05_duplication.tsv
