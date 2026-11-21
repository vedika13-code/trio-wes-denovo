#!/usr/bin/env bash
# Step 2: read quality control on the raw FASTQs
set -euo pipefail
source config.sh
mkdir -p $WORK/qc/fastqc_raw $RESULTS/qc

fastqc -t $THREADS -o $WORK/qc/fastqc_raw $DATA/fastq/*.fastq.gz
multiqc -f -o $RESULTS/qc -n 02_raw_multiqc $WORK/qc/fastqc_raw

# one line per file: PASS/WARN/FAIL for each FastQC module
echo -e "file\tmodule\tstatus" > $RESULTS/qc/02_fastqc_summary.tsv
for z in $WORK/qc/fastqc_raw/*_fastqc.zip; do
  unzip -p $z '*/summary.txt' | awk -F'\t' '{print $3"\t"$2"\t"$1}' >> $RESULTS/qc/02_fastqc_summary.tsv
done
grep -E "WARN|FAIL" $RESULTS/qc/02_fastqc_summary.tsv || echo "All modules PASS"
