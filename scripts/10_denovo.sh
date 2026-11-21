#!/usr/bin/env bash
# Step 10: find de novo candidates (child het, both parents hom-ref) with stepwise filters.
# Sample order in the VCF is fixed by step 9: [0]=child, [1]=father, [2]=mother.
set -euo pipefail
source config.sh
IN=$WORK/vcf/trio.final.vcf.gz
D=$WORK/denovo
mkdir -p $D $RESULTS/denovo
F=$RESULTS/denovo/10_funnel.tsv
echo -e "step\tdescription\tremaining" > $F
count() { bcftools view -H $1 | wc -l; }

echo -e "0\tall variant records\t$(count $IN)" >> $F

bcftools view -i 'GT[0]="het" && GT[1]="RR" && GT[2]="RR" && ALT!="*"' $IN -Oz -o $D/s1.vcf.gz
echo -e "1\tchild het, both parents hom-ref (v1 stopped here)\t$(count $D/s1.vcf.gz)" >> $F

bcftools view -f PASS $D/s1.vcf.gz -Oz -o $D/s2.vcf.gz
echo -e "2\tpassed hard filters\t$(count $D/s2.vcf.gz)" >> $F

bcftools view -i 'FMT/DP[0]>=10 && FMT/DP[1]>=10 && FMT/DP[2]>=10' $D/s2.vcf.gz -Oz -o $D/s3.vcf.gz
echo -e "3\tdepth >=10 in all three\t$(count $D/s3.vcf.gz)" >> $F

bcftools view -i 'FMT/GQ[0]>=20 && FMT/GQ[1]>=20 && FMT/GQ[2]>=20' $D/s3.vcf.gz -Oz -o $D/s4.vcf.gz
echo -e "4\tgenotype quality >=20 in all three\t$(count $D/s4.vcf.gz)" >> $F

bcftools view -i 'FMT/AD[0:1]/(FMT/AD[0:0]+FMT/AD[0:1])>=0.25 && FMT/AD[0:1]/(FMT/AD[0:0]+FMT/AD[0:1])<=0.75' \
  $D/s4.vcf.gz -Oz -o $D/s5.vcf.gz
echo -e "5\tchild allele balance 0.25-0.75\t$(count $D/s5.vcf.gz)" >> $F

bcftools view -i 'FMT/AD[1:1]=0 && FMT/AD[2:1]=0' $D/s5.vcf.gz -Oz -o $D/candidates.vcf.gz
tabix -f -p vcf $D/candidates.vcf.gz
echo -e "6\tzero alt reads in both parents\t$(count $D/candidates.vcf.gz)" >> $F

# readable table of the final candidates: genotype, read counts (ref,alt), depth, GQ for each person
echo -e "chrom\tpos\tref\talt\tchild_GT\tchild_AD\tchild_GQ\tfather_GT\tfather_AD\tfather_GQ\tmother_GT\tmother_AD\tmother_GQ" \
  > $RESULTS/denovo/10_candidates.tsv
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%AD\t%GQ]\n' $D/candidates.vcf.gz >> $RESULTS/denovo/10_candidates.tsv

column -t -s $'\t' $F 2>/dev/null || cat $F
echo; cat $RESULTS/denovo/10_candidates.tsv
