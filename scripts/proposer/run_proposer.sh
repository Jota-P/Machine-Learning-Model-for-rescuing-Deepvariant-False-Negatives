#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage:
  $0 --sample SAMPLE --bam BAM --ref REF --outdir OUTDIR \\
     --min-mapq N --min-bq N \\
     --min-snp-support N --min-indel-support N \\
     --max-indel-len N --max-depth N \\
     [--max-jobs N]

Required inputs:
  --sample SAMPLE              Sample name, e.g. HG007
  --bam BAM                    Input BAM file
  --ref REF                    Reference FASTA
  --outdir OUTDIR              Output directory

Required proposer parameters:
  --min-mapq N                 Minimum read mapping quality
  --min-bq N                   Minimum base quality
  --min-snp-support N          Minimum reads supporting a SNP candidate
  --min-indel-support N        Minimum reads supporting an INDEL candidate
  --max-indel-len N            Maximum insertion/deletion length to keep
  --max-depth N                Maximum mpileup depth

Optional:
  --max-jobs N                 Number of chromosomes processed in parallel [default: 4]
  -h, --help                   Show this help message
EOF
}

SAMPLE=""
BAM=""
REF=""
OUTBASE=""

MAX_JOBS=4
MIN_MAPQ=1
MIN_BQ=1
MIN_SNP_SUPPORT=1
MIN_INDEL_SUPPORT=1
MAX_INDEL_LEN=80
MAX_DEPTH=7000

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample) SAMPLE="$2"; shift 2 ;;
        --bam) BAM="$2"; shift 2 ;;
        --ref) REF="$2"; shift 2 ;;
        --outdir) OUTBASE="$2"; shift 2 ;;
        --max-jobs) MAX_JOBS="$2"; shift 2 ;;
        --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
        --min-bq) MIN_BQ="$2"; shift 2 ;;
        --min-snp-support) MIN_SNP_SUPPORT="$2"; shift 2 ;;
        --min-indel-support) MIN_INDEL_SUPPORT="$2"; shift 2 ;;
        --max-indel-len) MAX_INDEL_LEN="$2"; shift 2 ;;
        --max-depth) MAX_DEPTH="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument $1"; usage; exit 1 ;;
    esac
done

[[ -n "$SAMPLE" ]] || { echo "ERROR: --sample is required"; usage; exit 1; }
[[ -n "$BAM" ]] || { echo "ERROR: --bam is required"; usage; exit 1; }
[[ -n "$REF" ]] || { echo "ERROR: --ref is required"; usage; exit 1; }
[[ -n "$OUTBASE" ]] || { echo "ERROR: --outdir is required"; usage; exit 1; }

command -v samtools >/dev/null || { echo "ERROR: samtools not found"; exit 1; }
command -v python >/dev/null || { echo "ERROR: python not found"; exit 1; }

[[ -f "$BAM" ]] || { echo "ERROR: missing BAM: $BAM"; exit 1; }
[[ -f "$REF" ]] || { echo "ERROR: missing reference FASTA: $REF"; exit 1; }

LOGDIR="${OUTBASE}/logs"
mkdir -p "$OUTBASE" "$LOGDIR"

# Create indexes if missing
[[ -f "${BAM}.bai" ]] || samtools index -@ 4 "$BAM"
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

MANIFEST="${OUTBASE}/run_manifest.txt"

{
    echo "sample=${SAMPLE}"
    echo "bam=${BAM}"
    echo "ref=${REF}"
    echo "outdir=${OUTBASE}"
    echo "chromosomes=${CHROMS[*]}"
    echo "max_jobs=${MAX_JOBS}"
    echo "min_mapq=${MIN_MAPQ}"
    echo "min_bq=${MIN_BQ}"
    echo "min_snp_support=${MIN_SNP_SUPPORT}"
    echo "min_indel_support=${MIN_INDEL_SUPPORT}"
    echo "max_indel_len=${MAX_INDEL_LEN}"
    echo "max_depth=${MAX_DEPTH}"
    echo "samtools_version=$(samtools --version | head -n 1)"
    echo "python_version=$(python --version)"
    echo "run_date=$(date -Iseconds)"
} > "$MANIFEST"

run_chr() {
    local CHR="$1"

    local OUTVCF="${OUTBASE}/${SAMPLE}.chr${CHR}.candidates.vcf"
    local OUTBED="${OUTBASE}/${SAMPLE}.chr${CHR}.candidates.bed"
    local LOG="${LOGDIR}/${SAMPLE}.chr${CHR}.log"

    {
        echo "[$(date -Iseconds)] START chr${CHR}"

        if [[ -s "$OUTVCF" && -s "$OUTBED" ]]; then
            echo "[$(date -Iseconds)] SKIP chr${CHR}: outputs already exist"
            return 0
        fi

        python scripts/propose_candidates.py \
            --chrom "$CHR" \
            --bam "$BAM" \
            --ref "$REF" \
            --out-vcf "$OUTVCF" \
            --out-bed "$OUTBED" \
            --min-mapq "$MIN_MAPQ" \
            --min-bq "$MIN_BQ" \
            --min-snp-support "$MIN_SNP_SUPPORT" \
            --min-indel-support "$MIN_INDEL_SUPPORT" \
            --max-indel-len "$MAX_INDEL_LEN" \
            --max-depth "$MAX_DEPTH"

        echo "[$(date -Iseconds)] DONE chr${CHR}"
    } > "$LOG" 2>&1
}

for CHR in "${CHROMS[@]}"; do
    echo "Launching chr${CHR}"
    run_chr "$CHR" &

    while (( $(jobs -pr | wc -l) >= MAX_JOBS )); do
        wait -n
    done
done

wait

echo "All chromosomes finished."
echo "Manifest: $MANIFEST"
echo "Logs: $LOGDIR"
