#!/usr/bin/env bash
# Step 4: prepare the reference and align trimmed reads with BWA-MEM
set -euo pipefail
source config.sh
FA=$REF/chr$CHROM.fa
mkdir -p $WORK/bam $RESULTS/align

# --- reference indexes (each tool needs its own) ---
[ -f $FA.bwt ] || bwa index $FA                     # BWA's search index
[ -f $FA.fai ] || samtools faidx $FA                # samtools/GATK: fast lookup by position
[ -f ${FA%.fa}.dict ] || gatk CreateSequenceDictionary -R $FA -O ${FA%.fa}.dict   # GATK: chromosome names + lengths

# --- align each sample, sort, index ---
echo -e "sample\ttotal_reads\tmapped_pct\tproperly_paired_pct" > $RESULTS/align/04_flagstat.tsv
for S in $SAMPLES; do
  echo "== aligning $S =="
  RG="@RG\tID:$S\tSM:$S\tLB:${S}_lib\tPL:ILLUMINA\tPU:${S}_unit"
  bwa mem -t $THREADS -R "$RG" $FA \
      $WORK/trimmed/${S}_R1.paired.fq.gz $WORK/trimmed/${S}_R2.paired.fq.gz 2> $WORK/bam/$S.bwa.log \
    | samtools sort -@ $THREADS -o $WORK/bam/$S.sorted.bam -
  samtools index $WORK/bam/$S.sorted.bam
  samtools flagstat $WORK/bam/$S.sorted.bam > $RESULTS/align/$S.flagstat.txt
  total=$(awk 'NR==1{print $1}' $RESULTS/align/$S.flagstat.txt)
  mapped=$(grep -m1 "mapped (" $RESULTS/align/$S.flagstat.txt | sed 's/.*(\([0-9.]*%\).*/\1/')
  paired=$(grep "properly paired" $RESULTS/align/$S.flagstat.txt | sed 's/.*(\([0-9.]*%\).*/\1/')
  echo -e "$S\t$total\t$mapped\t$paired" >> $RESULTS/align/04_flagstat.tsv
done
cat $RESULTS/align/04_flagstat.tsv
