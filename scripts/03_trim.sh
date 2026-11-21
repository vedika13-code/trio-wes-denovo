#!/usr/bin/env bash
# Step 3: adapter + quality trimming with Trimmomatic, then QC again
set -euo pipefail
source config.sh
mkdir -p $WORK/trimmed $WORK/qc/fastqc_trimmed $RESULTS/qc

ADAPTERS=$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE-2.fa
echo -e "sample\tinput_pairs\tboth_surviving\tforward_only\treverse_only\tdropped" > $RESULTS/qc/03_trim_stats.tsv

for S in $SAMPLES; do
  echo "== trimming $S =="
  trimmomatic PE -threads $THREADS \
    $DATA/fastq/${S}_R1.fastq.gz $DATA/fastq/${S}_R2.fastq.gz \
    $WORK/trimmed/${S}_R1.paired.fq.gz $WORK/trimmed/${S}_R1.unpaired.fq.gz \
    $WORK/trimmed/${S}_R2.paired.fq.gz $WORK/trimmed/${S}_R2.unpaired.fq.gz \
    ILLUMINACLIP:$ADAPTERS:2:30:10:2:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50 \
    2> $WORK/trimmed/${S}.trimmomatic.log
  # pull the summary line out of the log into a table
  grep "Input Read Pairs" $WORK/trimmed/${S}.trimmomatic.log \
    | awk -v s=$S '{print s"\t"$4"\t"$7" "$8"\t"$12" "$13"\t"$17" "$18"\t"$20" "$21}' >> $RESULTS/qc/03_trim_stats.tsv
done

fastqc -t $THREADS -o $WORK/qc/fastqc_trimmed $WORK/trimmed/*.paired.fq.gz
multiqc -f -o $RESULTS/qc -n 03_trimmed_multiqc $WORK/qc/fastqc_trimmed $WORK/trimmed
cat $RESULTS/qc/03_trim_stats.tsv
