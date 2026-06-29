#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage:
  $0 --sample SAMPLE \\
     --labeled-tsv LABELED_TSV_GZ \\
     --bam BAM \\
     --ref REF_FASTA \\
     --outdir OUTDIR \\
     --min-mapq N \\
     --min-bq N \\
     --max-depth N \\
     [options]

Required:
  --sample SAMPLE              Sample/run name, e.g. HG007_100x
  --labeled-tsv FILE           Labeled candidate TSV file, compressed and indexed
  --bam FILE                   Input BAM file
  --ref FILE                   Reference FASTA file
  --outdir DIR                 Output directory for feature files
  --min-mapq N                 Minimum read mapping quality
  --min-bq N                   Minimum base quality
  --max-depth N                Maximum mpileup depth

Optional:
  --neg-keep-prob FLOAT        Probability of keeping negative examples [default: 1.0]
  --seed N                     Random seed for negative downsampling [default: 42]
  --max-jobs N                 Number of chromosomes processed in parallel [default: 1]
  --chroms "LIST"              Chromosomes to process [default: "1 2 ... 22"]
  --tmpdir DIR                 Temporary directory [default: system temp]
  -h, --help                   Show this help message

Output:
  One compressed feature TSV per chromosome.

Feature columns:
  CHROM POS REF ALT TYPE SUPPORT LABEL
  DP ALT_COUNT REF_COUNT VAF ALT_FWD ALT_REV SB ALT_BQ_MEAN
  CTX5 HOMOPOLY GC11
EOF
}

SAMPLE=""
LABELED=""
BAM=""
REF=""
OUTDIR=""

MIN_MAPQ=""
MIN_BQ=""
MAX_DEPTH=""

NEG_KEEP_PROB="1.0"
SEED="42"
MAX_JOBS="1"
CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
TMPDIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample) SAMPLE="$2"; shift 2 ;;
        --labeled-tsv) LABELED="$2"; shift 2 ;;
        --bam) BAM="$2"; shift 2 ;;
        --ref) REF="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
        --min-bq) MIN_BQ="$2"; shift 2 ;;
        --max-depth) MAX_DEPTH="$2"; shift 2 ;;
        --neg-keep-prob) NEG_KEEP_PROB="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        --max-jobs) MAX_JOBS="$2"; shift 2 ;;
        --chroms) CHROMS="$2"; shift 2 ;;
        --tmpdir) TMPDIR="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1"; usage; exit 1 ;;
    esac
done

[[ -n "$SAMPLE" ]] || { echo "ERROR: --sample is required"; usage; exit 1; }
[[ -n "$LABELED" ]] || { echo "ERROR: --labeled-tsv is required"; usage; exit 1; }
[[ -n "$BAM" ]] || { echo "ERROR: --bam is required"; usage; exit 1; }
[[ -n "$REF" ]] || { echo "ERROR: --ref is required"; usage; exit 1; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required"; usage; exit 1; }
[[ -n "$MIN_MAPQ" ]] || { echo "ERROR: --min-mapq is required"; usage; exit 1; }
[[ -n "$MIN_BQ" ]] || { echo "ERROR: --min-bq is required"; usage; exit 1; }
[[ -n "$MAX_DEPTH" ]] || { echo "ERROR: --max-depth is required"; usage; exit 1; }

command -v python >/dev/null || { echo "ERROR: python not found"; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not found"; exit 1; }
command -v tabix >/dev/null || { echo "ERROR: tabix not found"; exit 1; }

[[ -f "$LABELED" ]] || { echo "ERROR: missing labeled TSV: $LABELED"; exit 1; }
[[ -f "$BAM" ]] || { echo "ERROR: missing BAM: $BAM"; exit 1; }
[[ -f "$REF" ]] || { echo "ERROR: missing reference FASTA: $REF"; exit 1; }

mkdir -p "$OUTDIR"
LOGDIR="${OUTDIR}/logs"
mkdir -p "$LOGDIR"

if [[ -n "$TMPDIR" ]]; then
    mkdir -p "$TMPDIR"
    WORKROOT="$(mktemp -d "${TMPDIR}/features_${SAMPLE}.XXXXXX")"
else
    WORKROOT="$(mktemp -d)"
fi

cleanup() {
    rm -rf "$WORKROOT"
}
trap cleanup EXIT

# Index BAM if needed
if [[ ! -f "${BAM}.bai" ]]; then
    echo "Indexing BAM..."
    samtools index -@ 2 "$BAM"
fi

# Index reference if needed
if [[ ! -f "${REF}.fai" ]]; then
    echo "Indexing reference FASTA..."
    samtools faidx "$REF"
fi

# Index labeled TSV if needed
if [[ ! -f "${LABELED}.tbi" ]]; then
    echo "Indexing labeled TSV..."
    tabix -f -s 1 -b 2 -e 2 "$LABELED"
fi

MANIFEST="${OUTDIR}/feature_extraction_manifest.txt"

{
    echo "sample=${SAMPLE}"
    echo "labeled_tsv=${LABELED}"
    echo "bam=${BAM}"
    echo "ref=${REF}"
    echo "outdir=${OUTDIR}"
    echo "chroms=${CHROMS}"
    echo "min_mapq=${MIN_MAPQ}"
    echo "min_bq=${MIN_BQ}"
    echo "max_depth=${MAX_DEPTH}"
    echo "neg_keep_prob=${NEG_KEEP_PROB}"
    echo "seed=${SEED}"
    echo "max_jobs=${MAX_JOBS}"
    echo "samtools_version=$(samtools --version | head -n 1)"
    echo "python_version=$(python --version)"
    echo "run_date=$(date -Iseconds)"
} > "$MANIFEST"

run_chr() {
    local CHR="$1"
    local OUTFEAT="${OUTDIR}/${SAMPLE}.chr${CHR}.features.tsv.gz"
    local LOG="${LOGDIR}/${SAMPLE}.chr${CHR}.features.log"

    {
        echo "[$(date -Iseconds)] START chr${CHR}"
        echo "OUT=${OUTFEAT}"

        if [[ -s "$OUTFEAT" ]]; then
            echo "[$(date -Iseconds)] SKIP chr${CHR}: output already exists"
            return 0
        fi

        python scripts/extract_features.py \
            --chrom "$CHR" \
            --bam "$BAM" \
            --ref "$REF" \
            --labeled-tsv-gz "$LABELED" \
            --out-tsv-gz "$OUTFEAT" \
            --min-mapq "$MIN_MAPQ" \
            --min-bq "$MIN_BQ" \
            --max-depth "$MAX_DEPTH" \
            --neg-keep-prob "$NEG_KEEP_PROB" \
            --seed "$SEED" \
            --tmpdir "$WORKROOT"

        echo "[$(date -Iseconds)] DONE chr${CHR}"
    } > "$LOG" 2>&1
}

for CHR in $CHROMS; do
    echo "Launching chr${CHR}"
    run_chr "$CHR" &

    while (( $(jobs -pr | wc -l) >= MAX_JOBS )); do
        wait -n
    done
done

wait

echo
echo "Done. Feature files are in: $OUTDIR"
echo "Manifest: $MANIFEST"
echo "Logs: $LOGDIR"
ls -lh "$OUTDIR"/*.features.tsv.gz 2>/dev/null || true
