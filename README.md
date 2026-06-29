# Machine-Learning-Model-for-Rescuing-DeepVariant-False-Negatives

Machine learning pipeline for rescuing DeepVariant false-negative variant calls using candidate generation, candidate labelling, feature extraction, and XGBoost classification.

## Overview

This repository contains a variant rescue pipeline composed of four main steps:

1. **Candidate generation** from an aligned BAM file using `samtools mpileup`
2. **Candidate labelling** using a truth VCF and confident-region BED file
3. **Feature extraction** from pileup statistics and local sequence context
4. **XGBoost classification** to train SNP and INDEL rescue models

The goal is to generate a broad set of possible variant candidates and train machine learning models to identify candidates that correspond to DeepVariant false negatives.

## Repository structure

```text
scripts/
  ├── proposer/
    ├── run_proposer.sh
    ├── propose_candidates.py
    ├── label_candidates.sh
    ├── extract_features.sh
    ├── extract_features.py
  ├── xgboost
    └── train_xgboost.py
```

## Requirements

The pipeline requires:

* Python 3.10 or later
* samtools
* bcftools
* htslib
* pysam
* pandas
* numpy
* scikit-learn
* xgboost

A conda environment can be created with:

```bash
conda create -n dv-fn-rescue -c conda-forge -c bioconda \
  python=3.10 \
  samtools \
  bcftools \
  htslib \
  pysam \
  pandas \
  numpy \
  scikit-learn \
  xgboost
```

Then activate the environment:

```bash
conda activate dv-fn-rescue
```

Alternatively, the following `environment.yml` can be used:

```yaml
name: dv-fn-rescue
channels:
  - conda-forge
  - bioconda
dependencies:
  - python>=3.10
  - samtools
  - bcftools
  - htslib
  - pysam
  - pandas
  - numpy
  - scikit-learn
  - xgboost
```

Create the environment with:

```bash
conda env create -f environment.yml
conda activate dv-fn-rescue
```

## Input files

The pipeline requires:

* An aligned BAM file
* A reference FASTA file
* A truth VCF file
* A confident-region BED file

The BAM file should be indexed. If the BAM index is missing, the wrapper scripts can create it automatically using:

```bash
samtools index
```

The reference FASTA should also be indexed. If the FASTA index is missing, it can be created using:

```bash
samtools faidx
```

---

# 1. Candidate Variant Proposer

The candidate proposer generates permissive candidate variant sites from an aligned BAM file using `samtools mpileup`.

It scans each chromosome independently and writes candidate SNP and INDEL sites to VCF and BED files. The candidate proposer is intended to generate a broad set of possible variant sites that can later be filtered, labelled, or used as input features for a false-negative rescue model.

`run_proposer.sh` is the main wrapper script. It handles input validation, chromosome-level parallelization, logging, and run metadata.

`propose_candidates.py` contains the candidate-generation logic and parses the output of `samtools mpileup`.

## Usage

```bash
./scripts/run_proposer.sh \
  --sample SAMPLE_NAME \
  --bam /path/to/input.bam \
  --ref /path/to/reference.fa \
  --outdir results/SAMPLE_NAME/proposer \
  --min-mapq 10 \
  --min-bq 10 \
  --min-snp-support 2 \
  --min-indel-support 2 \
  --max-indel-len 50 \
  --max-depth 150 \
  --max-jobs 4
```

Replace the example threshold values with the values required for your experiment.

## Output

The proposer writes one VCF and one BED file per chromosome:

```text
results/SAMPLE_NAME/proposer/
├── SAMPLE_NAME.chr1.candidates.vcf
├── SAMPLE_NAME.chr1.candidates.bed
├── SAMPLE_NAME.chr2.candidates.vcf
├── SAMPLE_NAME.chr2.candidates.bed
├── ...
├── run_manifest.txt
└── logs/
```

---

# 2. Merge Candidate VCFs

The labelling step expects a single compressed and indexed candidate VCF.

The per-chromosome VCF files can be merged with `bcftools concat`:

```bash
bcftools concat \
  -Oz \
  -o results/SAMPLE_NAME/proposer/SAMPLE_NAME.candidates.merged.vcf.gz \
  results/SAMPLE_NAME/proposer/SAMPLE_NAME.chr*.candidates.vcf

tabix -f -p vcf results/SAMPLE_NAME/proposer/SAMPLE_NAME.candidates.merged.vcf.gz
```

---

# 3. Candidate Labelling

After candidate generation, candidate variants can be labelled using a truth VCF and a confident-region BED file.

Candidates are first restricted to confident regions. Then, each candidate is compared against the truth VCF using exact allele matching on:

```text
CHROM, POS, REF, ALT
```

Candidates that exactly match a truth variant receive label `1`. Candidates inside confident regions that do not match the truth VCF receive label `0`.

## Usage

```bash
./scripts/label_candidates.sh \
  --sample SAMPLE_NAME \
  --candidates-vcf results/SAMPLE_NAME/proposer/SAMPLE_NAME.candidates.merged.vcf.gz \
  --truth-vcf /path/to/truth.vcf.gz \
  --conf-bed /path/to/confident_regions.bed.gz \
  --out-tsv results/SAMPLE_NAME/SAMPLE_NAME.candidates.labeled.tsv.gz
```

## Output

```text
results/SAMPLE_NAME/
├── SAMPLE_NAME.candidates.labeled.tsv.gz
└── SAMPLE_NAME.candidates.labeled.tsv.gz.tbi
```

The labelled TSV contains:

| Column    | Description                                            |
| --------- | ------------------------------------------------------ |
| `CHROM`   | Chromosome                                             |
| `POS`     | 1-based variant position                               |
| `REF`     | Reference allele                                       |
| `ALT`     | Alternative allele                                     |
| `TYPE`    | Candidate type: SNP, INS, or DEL                       |
| `SUPPORT` | Alternative allele support from the candidate proposer |
| `LABEL`   | Truth label, where 1 is positive and 0 is negative     |

---

# 4. Feature Extraction

After candidate labelling, pileup-based features can be extracted for each labelled candidate site.

The feature extractor uses the labelled TSV, the original BAM file, and the reference FASTA. For each candidate site, it runs `samtools mpileup` and computes read-support, allele-balance, strand-bias, base-quality, and local sequence-context features.

## Usage

```bash
./scripts/extract_features.sh \
  --sample SAMPLE_NAME \
  --labeled-tsv results/SAMPLE_NAME/SAMPLE_NAME.candidates.labeled.tsv.gz \
  --bam /path/to/input.bam \
  --ref /path/to/reference.fa \
  --outdir results/SAMPLE_NAME/features \
  --min-mapq 10 \
  --min-bq 10 \
  --max-depth 5000 \
  --neg-keep-prob 1.0 \
  --seed 42 \
  --max-jobs 4
```

## Output

The script writes one compressed feature file per chromosome:

```text
results/SAMPLE_NAME/features/
├── SAMPLE_NAME.chr1.features.tsv.gz
├── SAMPLE_NAME.chr2.features.tsv.gz
├── ...
├── feature_extraction_manifest.txt
└── logs/
```

The feature files contain:

| Column        | Description                                       |
| ------------- | ------------------------------------------------- |
| `CHROM`       | Chromosome                                        |
| `POS`         | 1-based candidate position                        |
| `REF`         | Reference allele                                  |
| `ALT`         | Alternative allele                                |
| `TYPE`        | Candidate type: SNP, INS, or DEL                  |
| `SUPPORT`     | Support value from the candidate proposer         |
| `LABEL`       | Truth label                                       |
| `DP`          | Read depth at the candidate site                  |
| `ALT_COUNT`   | Number of reads supporting the alternative allele |
| `REF_COUNT`   | Number of reads supporting the reference allele   |
| `VAF`         | Alternative allele fraction                       |
| `ALT_FWD`     | Forward-strand alternative support                |
| `ALT_REV`     | Reverse-strand alternative support                |
| `SB`          | Strand-bias proxy                                 |
| `ALT_BQ_MEAN` | Mean base quality of SNP alternative observations |
| `CTX5`        | 5 bp local sequence context                       |
| `HOMOPOLY`    | Reference homopolymer length                      |
| `GC11`        | GC fraction in an 11 bp window                    |

---

# 5. Merge Feature Files

The XGBoost training script expects one feature file per sample. If feature extraction was run chromosome by chromosome, the chromosome files should be merged first.

Example:

```bash
zcat results/SAMPLE_NAME/features/SAMPLE_NAME.chr1.features.tsv.gz | head -n 1 \
  | gzip > results/SAMPLE_NAME/features/SAMPLE_NAME.rescue_features.allchr.tsv.gz

for chr in {1..22}; do
  zcat results/SAMPLE_NAME/features/SAMPLE_NAME.chr${chr}.features.tsv.gz | tail -n +2
done | gzip >> results/SAMPLE_NAME/features/SAMPLE_NAME.rescue_features.allchr.tsv.gz
```

---

# 6. Train XGBoost Rescue Classifiers

After feature extraction, SNP and INDEL rescue classifiers can be trained with XGBoost.

The training script expects compressed feature TSV files. Each input file can correspond to one sample or to a merged all-chromosome feature file for one sample.

The script trains two binary classifiers:

| Model                         | Variant types |
| ----------------------------- | ------------- |
| `xgb_snp_rescue_classifier`   | SNP           |
| `xgb_indel_rescue_classifier` | INS and DEL   |

## Usage

```bash
python scripts/train_xgboost.py \
  --train HG001=results/HG001/features/HG001.rescue_features.allchr.tsv.gz \
  --train HG002=results/HG002/features/HG002.rescue_features.allchr.tsv.gz \
  --train HG003=results/HG003/features/HG003.rescue_features.allchr.tsv.gz \
  --train HG004=results/HG004/features/HG004.rescue_features.allchr.tsv.gz \
  --train HG005=results/HG005/features/HG005.rescue_features.allchr.tsv.gz \
  --valid HG006=results/HG006/features/HG006.rescue_features.allchr.tsv.gz \
  --test HG007=results/HG007/features/HG007.rescue_features.allchr.tsv.gz \
  --outdir results/xgboost/classifier \
  --seed 13 \
  --nthread 6 \
  --models both
```

## Output

```text
results/xgboost/classifier/
├── training_manifest.json
├── xgb_snp_rescue_classifier.json
├── xgb_snp_rescue_classifier.metrics.json
├── xgb_snp_rescue_classifier.feature_importance.tsv
├── xgb_snp_rescue_classifier.valid_thresholds.tsv
├── xgb_snp_rescue_classifier.test_at_valid_selected_thresholds.tsv
├── xgb_snp_rescue_classifier.HG006.valid_predictions.tsv.gz
├── xgb_snp_rescue_classifier.HG007.test_predictions.tsv.gz
├── xgb_indel_rescue_classifier.json
├── xgb_indel_rescue_classifier.metrics.json
├── xgb_indel_rescue_classifier.feature_importance.tsv
├── xgb_indel_rescue_classifier.valid_thresholds.tsv
├── xgb_indel_rescue_classifier.test_at_valid_selected_thresholds.tsv
├── xgb_indel_rescue_classifier.HG006.valid_predictions.tsv.gz
└── xgb_indel_rescue_classifier.HG007.test_predictions.tsv.gz
```

---

# Notes

* Chromosomes are processed independently.
* The default chromosome list is `1` to `22`.
* Chromosome names must match the reference FASTA and BAM file.

  * GRCh37/hs37d5 usually uses `1`, `2`, ..., `22`.
  * References with UCSC-style names may use `chr1`, `chr2`, ..., `chr22`.
* The candidate proposer is intentionally permissive.
* The exact parameters used for each run are saved in manifest files.
* Large data files such as BAM, VCF, BED, and output TSV files should not be committed to GitHub.
