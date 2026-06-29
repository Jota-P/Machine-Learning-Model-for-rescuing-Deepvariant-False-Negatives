#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage:
  $0 --sample SAMPLE \\
     --candidates-vcf CANDIDATES_VCF_GZ \\
     --truth-vcf TRUTH_VCF_GZ \\
     --conf-bed CONFIDENT_BED_GZ \\
     --out-tsv OUTPUT_TSV_GZ

Required:
  --sample SAMPLE              Sample/run name, e.g. HG007_100x
  --candidates-vcf FILE        Candidate VCF file, compressed with bgzip
  --truth-vcf FILE             Truth VCF file, compressed with bgzip
  --conf-bed FILE              Confident-region BED file, optionally compressed
  --out-tsv FILE               Output labeled TSV file ending in .tsv.gz

Optional:
  --tmpdir DIR                 Temporary directory [default: system temp]
  -h, --help                   Show this help message

Output columns:
  CHROM POS REF ALT TYPE SUPPORT LABEL

LABEL:
  1 = candidate matches a truth variant exactly by CHROM, POS, REF, ALT
  0 = candidate is inside the confident region but does not match truth
EOF
}

SAMPLE=""
CAND_VCF=""
TRUTH_VCF=""
CONF_BED=""
OUT_TSV=""
TMPDIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample) SAMPLE="$2"; shift 2 ;;
        --candidates-vcf) CAND_VCF="$2"; shift 2 ;;
        --truth-vcf) TRUTH_VCF="$2"; shift 2 ;;
        --conf-bed) CONF_BED="$2"; shift 2 ;;
        --out-tsv) OUT_TSV="$2"; shift 2 ;;
        --tmpdir) TMPDIR="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1"; usage; exit 1 ;;
    esac
done

[[ -n "$SAMPLE" ]] || { echo "ERROR: --sample is required"; usage; exit 1; }
[[ -n "$CAND_VCF" ]] || { echo "ERROR: --candidates-vcf is required"; usage; exit 1; }
[[ -n "$TRUTH_VCF" ]] || { echo "ERROR: --truth-vcf is required"; usage; exit 1; }
[[ -n "$CONF_BED" ]] || { echo "ERROR: --conf-bed is required"; usage; exit 1; }
[[ -n "$OUT_TSV" ]] || { echo "ERROR: --out-tsv is required"; usage; exit 1; }

command -v bcftools >/dev/null || { echo "ERROR: bcftools not found"; exit 1; }
command -v tabix >/dev/null || { echo "ERROR: tabix not found"; exit 1; }
command -v bgzip >/dev/null || { echo "ERROR: bgzip not found"; exit 1; }

[[ -f "$CAND_VCF" ]] || { echo "ERROR: missing candidate VCF: $CAND_VCF"; exit 1; }
[[ -f "$TRUTH_VCF" ]] || { echo "ERROR: missing truth VCF: $TRUTH_VCF"; exit 1; }
[[ -f "$CONF_BED" ]] || { echo "ERROR: missing confident BED: $CONF_BED"; exit 1; }

OUTDIR="$(dirname "$OUT_TSV")"
mkdir -p "$OUTDIR"

if [[ -n "$TMPDIR" ]]; then
    mkdir -p "$TMPDIR"
    WORK="$(mktemp -d "${TMPDIR}/label_${SAMPLE}.XXXXXX")"
else
    WORK="$(mktemp -d)"
fi

cleanup() {
    rm -rf "$WORK"
}
trap cleanup EXIT

TMP_OUT="${OUT_TSV}.tmp"

echo "SAMPLE=${SAMPLE}"
echo "CANDIDATES_VCF=${CAND_VCF}"
echo "TRUTH_VCF=${TRUTH_VCF}"
echo "CONF_BED=${CONF_BED}"
echo "OUT_TSV=${OUT_TSV}"
echo "WORKDIR=${WORK}"
echo "bcftools=$(bcftools --version | head -n 1)"
echo "tabix=$(command -v tabix)"
echo "bgzip=$(command -v bgzip)"

# Make sure input VCFs are indexed.
if [[ ! -f "${CAND_VCF}.tbi" && ! -f "${CAND_VCF}.csi" ]]; then
    echo "Indexing candidate VCF..."
    tabix -f -p vcf "$CAND_VCF"
fi

if [[ ! -f "${TRUTH_VCF}.tbi" && ! -f "${TRUTH_VCF}.csi" ]]; then
    echo "Indexing truth VCF..."
    tabix -f -p vcf "$TRUTH_VCF"
fi

echo "Restricting candidates to confident regions..."

bcftools view \
    -R "$CONF_BED" \
    -Oz \
    -o "${WORK}/cand.conf.vcf.gz" \
    "$CAND_VCF"

tabix -f -p vcf "${WORK}/cand.conf.vcf.gz"

echo "Intersecting candidates with truth variants..."

bcftools isec \
    -c all \
    -p "${WORK}/isec" \
    "${WORK}/cand.conf.vcf.gz" \
    "$TRUTH_VCF" \
    >/dev/null

echo "Building truth-matching key set..."

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    "${WORK}/isec/0002.vcf" \
    | awk 'BEGIN{OFS="\t"} {print $1":"$2":"$3":"$4}' \
    > "${WORK}/truth_keys.txt"

echo "Writing labeled TSV..."

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/TYPE\t%INFO/SUPPORT\n' \
    "${WORK}/cand.conf.vcf.gz" \
    | awk '
        BEGIN {OFS="\t"}
        NR==FNR {
            truth[$1]=1
            next
        }
        {
            key=$1":"$2":"$3":"$4
            label=(key in truth) ? 1 : 0
            print $0, label
        }
    ' "${WORK}/truth_keys.txt" - \
    | bgzip -c > "$TMP_OUT"

tabix -f -s 1 -b 2 -e 2 "$TMP_OUT"

mv "$TMP_OUT" "$OUT_TSV"
mv "${TMP_OUT}.tbi" "${OUT_TSV}.tbi"

echo "DONE"
ls -lh "$OUT_TSV" "${OUT_TSV}.tbi"

echo -n "Rows: "
zcat "$OUT_TSV" | wc -l

echo -n "Positives: "
zcat "$OUT_TSV" | awk '$7==1{c++} END{print c+0}'

echo -n "Negatives: "
zcat "$OUT_TSV" | awk '$7==0{c++} END{print c+0}'
