#!/usr/bin/env bash
# Step 1: download chromosome-20 exome reads for the trio, the reference, and capture targets
set -euo pipefail          # stop at the first error instead of carrying on with bad files
source config.sh
mkdir -p $DATA/fastq $DATA/bam_slices $REF $RESULTS

GIAB=https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio
OSLO=OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample
declare -A URL=(
  [HG002]=$GIAB/HG002_NA24385_son/${OSLO}_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup
  [HG003]=$GIAB/HG003_NA24149_father/${OSLO}_HG003-EEogPU_v02-KIT-Av5_TCTTCACA_L008.posiSrt.markDup
  [HG004]=$GIAB/HG004_NA24143_mother/${OSLO}_HG004-EEogPU_v02-KIT-Av5_CCGAAGTA_L008.posiSrt.markDup
)

# 1. Reads: take only chromosome 20 from each remote BAM, then turn them back into FASTQ
for S in $SAMPLES; do
  echo "== $S: extracting chromosome $CHROM =="
  samtools view -b -X -o $DATA/bam_slices/$S.chr$CHROM.bam ${URL[$S]}.bam ${URL[$S]}.bai $CHROM
  samtools collate -u -O -T $DATA/bam_slices/${S}_tmp $DATA/bam_slices/$S.chr$CHROM.bam \
    | samtools fastq -n -1 $DATA/fastq/${S}_R1.fastq.gz -2 $DATA/fastq/${S}_R2.fastq.gz \
                        -0 /dev/null -s /dev/null -
done

# 2. Reference: chromosome 20 of GRCh37 (UCSC hg19 chr20 is the same sequence; rename to "20")
curl -fL https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr20.fa.gz \
  | gunzip -c | sed 's/^>chr20/>20/' > $REF/chr$CHROM.fa

# 3. Exome capture targets (Agilent SureSelect v5, GRCh37), chromosome 20 only
curl -fL https://storage.googleapis.com/deepvariant/exome-case-study-testdata/agilent_sureselect_human_all_exon_v5_b37_targets.bed \
  | awk -v c=$CHROM '$1==c' > $REF/targets.chr$CHROM.bed

# 4. Record how many read pairs each person has
echo -e "sample\tread_pairs" > $RESULTS/01_read_counts.tsv
for S in $SAMPLES; do
  n=$(( $(zcat $DATA/fastq/${S}_R1.fastq.gz | wc -l) / 4 ))
  echo -e "$S\t$n" >> $RESULTS/01_read_counts.tsv
done
cat $RESULTS/01_read_counts.tsv
echo "Targets on chr$CHROM: $(wc -l < $REF/targets.chr$CHROM.bed)"
