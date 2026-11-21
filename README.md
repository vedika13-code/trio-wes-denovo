# Trio exome de novo variant calling, validated against GIAB truth sets

A reproducible pipeline that calls **de novo variants** (present in a child, absent in both
parents) from trio whole-exome sequencing, and then **checks every candidate against
independent truth data**.

**Data:** the Genome in a Bottle (GIAB) Ashkenazim trio: HG002 (son), HG003 (father),
HG004 (mother). Exome sequencing from Oslo University Hospital (Agilent SureSelect v5),
2×126 bp Illumina reads. To keep it laptop-sized, the pipeline uses **chromosome 20** at
full depth. Reference: GRCh37.

## Results

**Coverage on chr20 exome targets** (duplicates and MAPQ < 20 removed)

| | HG002 | HG003 | HG004 |
|---|---|---|---|
| Mean depth | 199× | 169× | 186× |
| Target bases ≥ 20× | 99.1% | 98.8% | 99.1% |
| PCR duplicates | 4.2% | 3.6% | 4.2% |

**De novo filter funnel** (`results/denovo/10_funnel.tsv`)

| Step | Candidates |
|---|---|
| Child het, both parents hom-ref | 22 |
| + passed GATK hard filters | 15 |
| + depth ≥ 10 in all three | 9 |
| + genotype quality ≥ 20 in all three | 6 |
| + child allele balance 0.25–0.75 | 3 |
| + zero alt reads in both parents | 3 |

**Validation against GIAB v4.2.1 truth sets for all three trio members** (`results/denovo/11_*`)

| | Before filters (22) | After filters (3) |
|---|---|---|
| False call in child | 11 | 2 |
| Inherited (missed in a parent) | 2 | 0 |
| Outside GIAB confident regions | 9 | 1 |
| Confirmed de novo | 0 | 0 |

**Interpretation.** Roughly 1–2 de novo coding variants are expected per whole exome; chr20 is
about 2–3% of the exome, so the expected count here is close to zero. The filters removed both
inherited variants and 9 of 11 false child calls. The remaining false positives are **small
indels with low local depth**, the typical failure mode of short-read calling in repetitive
sequence.

## Pipeline

| Script | Step | Tools |
|---|---|---|
| `01_get_data.sh` | Download chr20 exome reads (remote BAM slice → FASTQ), reference, capture targets | samtools |
| `01b_get_truth.sh` | GIAB v4.2.1 truth VCFs + confident regions for all three | bcftools |
| `02_qc.sh` | Raw read QC | FastQC, MultiQC |
| `03_trim.sh` | Adapter and quality trimming | Trimmomatic |
| `04_align.sh` | Alignment with read groups, sort, index | BWA-MEM, samtools |
| `05_markdup.sh` | Mark PCR/optical duplicates | samtools markdup |
| `06_coverage.sh` | Depth over exome targets | samtools depth |
| `07_call_gvcf.sh` | Per-sample calling, GVCF mode, targets ± 100 bp | GATK HaplotypeCaller |
| `08_joint_genotype.sh` | Joint genotyping of the trio | GATK CombineGVCFs, GenotypeGVCFs |
| `09_filter.sh` | Hard filters (SNP/indel specific), normalisation | GATK, bcftools norm |
| `10_denovo.sh` | Stepwise de novo filters with counts | bcftools |
| `11_check_truth.sh` | Check candidates against GIAB truth for child, father, mother | bcftools |

## How to run

```bash
micromamba create -f environment.yml     # or: conda env create -f environment.yml
micromamba activate trio
./run_all.sh
```

Settings (samples, chromosome, folders, threads) are in `config.sh`. Large files (`data/`,
`ref/`, `work/`) are not tracked in git; `01_get_data.sh` recreates them.

## Design choices and limitations

- **One chromosome.** Reads were extracted from chr20 of the original alignments and realigned
  to chr20 only. A few reads from repeats with copies elsewhere may be misplaced; a
  whole-genome run would avoid this.
- **Hard filters, not VQSR.** VQSR needs many samples; GATK recommends hard filters for a
  single trio. Most filtered sites failed SOR (strand bias), which is common at exome target edges.
- **No BQSR** in this version.
- **samtools markdup** instead of GATK MarkDuplicates, which hit intermittent write failures in
  the development VM; both mark duplicates using the positions of both read mates.
- **Indels in repeats** are the main remaining false-positive source. A further filter (for
  example on local depth or repeat context) should be tested on a different chromosome or
  sample, not tuned on this truth set.

## Version history

v1 (`Trio-WES-RareDisease-Analysis`) used a 4M-read-pair subsample at about 0.5× depth and
reported 12,776 candidate de novo variants, mostly sites where parents had no reads but were
still genotyped 0/0. This version fixes that with adequate depth, quality filters in all three
samples, and validation against truth data.

## Software

samtools / bcftools / htslib 1.20 · BWA 0.7.18 · GATK 4.6.1.0 · Trimmomatic 0.39 ·
FastQC 0.12.1 · MultiQC 1.25

## Data sources

- GIAB Ashkenazim trio exome BAMs (Oslo University Hospital), NCBI GIAB FTP
- GIAB v4.2.1 small-variant benchmarks for HG002, HG003, HG004 (GRCh37)
- Agilent SureSelect v5 targets (b37), from the DeepVariant exome case study
- UCSC hg19 chr20 (sequence-identical to GRCh37 chromosome 20)

Author: Vedika Goyal, VIT Vellore
