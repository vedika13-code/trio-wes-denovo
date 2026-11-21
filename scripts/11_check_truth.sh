#!/usr/bin/env bash
# Step 11: check de novo candidates against GIAB v4.2.1 truth sets for child, father and mother
set -euo pipefail
source config.sh
FA="$REF/chr$CHROM.fa"
T="$DATA/truth"
mkdir -p "$RESULTS/denovo"

# normalise truth sets the same way as our calls, so indels are written identically
for S in $SAMPLES; do
  [ -s "$T/$S.norm.vcf.gz.tbi" ] && continue
  bcftools norm -m -any -f "$FA" "$T/$S.chr$CHROM.vcf.gz" -Oz -o "$T/$S.norm.vcf.gz" 2>/dev/null
  tabix -f -p vcf "$T/$S.norm.vcf.gz"
done

# is position $2 on chrom $1 inside sample $3's GIAB confident regions?
in_bed() {
  awk -v c="$1" -v p="$2" '$1==c && $2<p && $3>=p {f=1} END{print f?1:0}' "$T/$3.chr$CHROM.bed"
}
# does sample $5's truth set contain exactly this variant (chrom $1, pos $2, ref $3, alt $4)?
in_truth() {
  bcftools view -H -r "$1:$2-$2" "$T/$5.norm.vcf.gz" \
    | awk -v p="$2" -v r="$3" -v a="$4" '$2==p && $4==r && $5==a {f=1} END{print f?1:0}'
}

check() {   # $1 = VCF to check, $2 = output table
  echo -e "chrom\tpos\tref\talt\tchild_in_bed\tfather_in_bed\tmother_in_bed\tchild_has\tfather_has\tmother_has\tverdict" > "$2"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$1" | while read -r c p r a; do
    cb=$(in_bed "$c" "$p" "$CHILD");  fb=$(in_bed "$c" "$p" "$FATHER");  mb=$(in_bed "$c" "$p" "$MOTHER")
    ch=$(in_truth "$c" "$p" "$r" "$a" "$CHILD")
    fh=$(in_truth "$c" "$p" "$r" "$a" "$FATHER")
    mh=$(in_truth "$c" "$p" "$r" "$a" "$MOTHER")
    if   [ "$fh" = 1 ] || [ "$mh" = 1 ]; then v=inherited
    elif [ "$cb" = 1 ] && [ "$ch" = 0 ]; then v=false_call_in_child
    elif [ "$ch" = 1 ] && [ "$cb" = 1 ] && [ "$fb" = 1 ] && [ "$mb" = 1 ]; then v=TRUE_DE_NOVO
    else v=outside_GIAB_confident_regions; fi
    echo -e "$c\t$p\t$r\t$a\t$cb\t$fb\t$mb\t$ch\t$fh\t$mh\t$v" >> "$2"
  done
}

check "$WORK/denovo/s1.vcf.gz"         "$RESULTS/denovo/11_truth_before_filters.tsv"
check "$WORK/denovo/candidates.vcf.gz" "$RESULTS/denovo/11_truth_final.tsv"

echo "== Before quality filters (step 1: $(tail -n +2 "$RESULTS/denovo/11_truth_before_filters.tsv" | wc -l) sites) =="
cut -f11 "$RESULTS/denovo/11_truth_before_filters.tsv" | tail -n +2 | sort | uniq -c
echo "== Final candidates =="
cat "$RESULTS/denovo/11_truth_final.tsv"
