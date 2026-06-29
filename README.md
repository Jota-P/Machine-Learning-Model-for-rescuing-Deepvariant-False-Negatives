# Machine-Learning-Model-for-rescuing-Deepvariant-False-Negatives
Machine learning pipeline for rescuing DeepVariant false-negative variant calls using candidate generation, feature extraction, and XGBoost classification.

# Candidate Variant Proposer

This script generates permissive candidate variant sites from an aligned BAM file using `samtools mpileup`.

It scans each chromosome independently and writes candidate SNP and INDEL sites to VCF and BED files. The candidate proposer is intended to generate a broad set of possible variant sites that can later be filtered, labelled, or used as input features for a false-negative rescue model.

run_proposer.sh is the main wrapper script. It handles input validation, chromosome-level parallelization, logging, and run metadata.
propose_candidates.py contains the candidate-generation logic and parses the output of samtools mpileup.


# Requirements
  Python 3.10 or later
  samtools

A conda environment can be created with:

```bash
conda create -n dv-fn-rescue -c conda-forge -c bioconda python=3.10 samtools
conda activate dv-fn-rescue
```

The following required input files are needed:
  - An aligned BAM file
  - A reference FASTA file

The BAM file should be indexed. If the BAM index is missing, the wrapper script will create it automatically using:
```bash
samtools index
```

The reference FASTA should also be indexed. If the FASTA index is missing, the wrapper script will create it automatically using:
  ```bash
samtools faidx
```


#### Usage

```bash
./scripts/run_proposer.sh \
  --sample SAMPLE_NAME \
  --bam /path/to/input.bam \
  --ref /path/to/reference.fa \
  --outdir /path/to/output_directory \
  --min-mapq {var} \
  --min-bq {var} \
  --min-snp-support {var} \
  --min-indel-support {var} \
  --max-indel-len {var} \
  --max-depth {var} \
  --max-jobs {var}
```

# Labelling

After candidate generation, candidate variants can be labelled using a truth VCF and a confident-region BED file.

Candidates are first restricted to the confident regions. Then, each candidate is compared against the truth VCF using exact allele matching.
This labelling expects a merged VCF from the proposer outputs. Which can be done via:

```bash
bcftools concat
```

#### Usage
```bash
./scripts/label_candidates.sh \
  --sample {var} \
  --candidates-vcf {var} \
  --truth-vcf /path/to/{TRUTH}.vcf.gz \
  --conf-bed /path/to/{CONF}.bed.gz \
  --out-tsv results/{SAMPLE}/{OUT}.labeled.tsv.gz
```




