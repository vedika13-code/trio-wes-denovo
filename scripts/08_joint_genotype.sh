#!/usr/bin/env bash
# Step 8: combine the three GVCFs and genotype the trio together
set -euo pipefail
source config.sh
FA=$REF/chr$CHROM.fa
BED=$REF/targets.chr$CHROM.bed
mkdir -p $WORK/vcf $RESULTS/vcf

V=""; for S in $SAMPLES; do V="$V -V $WORK/gvcf/$S.g.vcf.gz"; done

gatk --java-options "-Xmx6g" CombineGVCFs -R $FA $V -L $BED -ip 100 \
  --tmp-dir $WORK/tmp -O $WORK/vcf/trio.combined.g.vcf.gz 2> $WORK/vcf/combine.log

gatk --java-options "-Xmx6g" GenotypeGVCFs -R $FA -V $WORK/vcf/trio.combined.g.vcf.gz -L $BED -ip 100 \
  --tmp-dir $WORK/tmp -O $WORK/vcf/trio.raw.vcf.gz 2> $WORK/vcf/genotype.log

# summary: number of variant sites, split into SNPs and indels
{
  echo -e "metric\tcount"
  echo -e "all_sites\t$(bcftools view -H $WORK/vcf/trio.raw.vcf.gz | wc -l)"
  echo -e "snps\t$(bcftools view -H -v snps $WORK/vcf/trio.raw.vcf.gz | wc -l)"
  echo -e "indels\t$(bcftools view -H -v indels $WORK/vcf/trio.raw.vcf.gz | wc -l)"
} > $RESULTS/vcf/08_joint_raw_counts.tsv
cat $RESULTS/vcf/08_joint_raw_counts.tsv
