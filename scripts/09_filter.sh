#!/usr/bin/env bash
# Step 9: GATK hard filters (separate thresholds for SNPs and indels), then normalise
set -euo pipefail
source config.sh
FA=$REF/chr$CHROM.fa
IN=$WORK/vcf/trio.raw.vcf.gz
mkdir -p $RESULTS/vcf

# SNPs
gatk SelectVariants -V $IN -select-type SNP -O $WORK/vcf/snps.raw.vcf.gz 2> /dev/null
gatk VariantFiltration -V $WORK/vcf/snps.raw.vcf.gz -O $WORK/vcf/snps.filt.vcf.gz \
  -filter "QD < 2.0"              --filter-name "QD2" \
  -filter "QUAL < 30.0"           --filter-name "QUAL30" \
  -filter "SOR > 3.0"             --filter-name "SOR3" \
  -filter "FS > 60.0"             --filter-name "FS60" \
  -filter "MQ < 40.0"             --filter-name "MQ40" \
  -filter "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" 2> /dev/null

# Indels (and sites mixing SNP + indel alleles)
gatk SelectVariants -V $IN -select-type INDEL -select-type MIXED -O $WORK/vcf/indels.raw.vcf.gz 2> /dev/null
gatk VariantFiltration -V $WORK/vcf/indels.raw.vcf.gz -O $WORK/vcf/indels.filt.vcf.gz \
  -filter "QD < 2.0"               --filter-name "QD2" \
  -filter "QUAL < 30.0"            --filter-name "QUAL30" \
  -filter "FS > 200.0"             --filter-name "FS200" \
  -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" 2> /dev/null

gatk MergeVcfs -I $WORK/vcf/snps.filt.vcf.gz -I $WORK/vcf/indels.filt.vcf.gz \
  -O $WORK/vcf/trio.filtered.vcf.gz 2> /dev/null

# Normalise: fixed sample order (child, father, mother), one alt allele per line, left-aligned indels
bcftools view -s $CHILD,$FATHER,$MOTHER $WORK/vcf/trio.filtered.vcf.gz \
  | bcftools norm -m -any -f $FA -Oz -o $WORK/vcf/trio.final.vcf.gz
tabix -f -p vcf $WORK/vcf/trio.final.vcf.gz

{
  echo -e "metric\tcount"
  echo -e "records_after_split\t$(bcftools view -H $WORK/vcf/trio.final.vcf.gz | wc -l)"
  echo -e "PASS\t$(bcftools view -H -f PASS $WORK/vcf/trio.final.vcf.gz | wc -l)"
  echo -e "filtered_out\t$(bcftools view -H -e 'FILTER="PASS"' $WORK/vcf/trio.final.vcf.gz | wc -l)"
} > $RESULTS/vcf/09_filter_counts.tsv
echo "-- which filters failed (count per filter name) --" >> $RESULTS/vcf/09_filter_counts.tsv
bcftools query -i 'FILTER!="PASS"' -f '%FILTER\n' $WORK/vcf/trio.final.vcf.gz | tr ';' '\n' | sort | uniq -c >> $RESULTS/vcf/09_filter_counts.tsv
cat $RESULTS/vcf/09_filter_counts.tsv
