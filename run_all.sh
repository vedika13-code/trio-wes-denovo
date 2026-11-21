#!/usr/bin/env bash
# Run the full trio de novo pipeline, step by step. Each step stops on error.
set -euo pipefail
for s in 01_get_data 01b_get_truth 02_qc 03_trim 04_align 05_markdup \
         06_coverage 07_call_gvcf 08_joint_genotype 09_filter 10_denovo 11_check_truth; do
  echo "########## $s ##########"
  bash scripts/$s.sh
done
echo "Done. Key results: results/denovo/10_funnel.tsv and results/denovo/11_truth_final.tsv"
