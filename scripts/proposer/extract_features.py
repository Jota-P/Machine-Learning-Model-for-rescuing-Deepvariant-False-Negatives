#!/usr/bin/env python3

import argparse
import gzip
import os
import random
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pysam


HEADER = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "TYPE",
    "SUPPORT",
    "LABEL",
    "DP",
    "ALT_COUNT",
    "REF_COUNT",
    "VAF",
    "ALT_FWD",
    "ALT_REV",
    "SB",
    "ALT_BQ_MEAN",
    "CTX5",
    "HOMOPOLY",
    "GC11",
]

INDEL_RE = re.compile(r"([+-])(\d+)")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract pileup-based features for labeled candidate variants."
    )

    parser.add_argument("--chrom", required=True, help="Chromosome/contig to process")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--ref", required=True, help="Reference FASTA")
    parser.add_argument(
        "--labeled-tsv-gz",
        required=True,
        help="Compressed and tabix-indexed labeled candidate TSV",
    )
    parser.add_argument(
        "--out-tsv-gz",
        required=True,
        help="Output compressed feature TSV",
    )

    parser.add_argument("--min-mapq", type=int, required=True)
    parser.add_argument("--min-bq", type=int, required=True)
    parser.add_argument("--max-depth", type=int, required=True)
    parser.add_argument("--neg-keep-prob", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--tmpdir", default=None)

    return parser.parse_args()


def write_header_only(path):
    with gzip.open(path, "wt") as out:
        out.write("\t".join(HEADER) + "\n")


def read_labeled_rows(args):
    cmd = ["tabix", args.labeled_tsv_gz, args.chrom]

    process = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if process.returncode != 0 and process.stdout.strip() == "":
        return []

    rows = []

    for line in process.stdout.splitlines():
        if not line.strip():
            continue

        parts = line.rstrip("\n").split("\t")

        if len(parts) < 7:
            continue

        rows.append(
            {
                "CHROM": parts[0],
                "POS": int(parts[1]),
                "REF": parts[2],
                "ALT": parts[3],
                "TYPE": parts[4],
                "SUPPORT": int(parts[5]),
                "LABEL": int(parts[6]),
            }
        )

    return rows


def downsample_negatives(rows, neg_keep_prob, seed, chrom):
    if neg_keep_prob >= 1.0:
        return rows

    if neg_keep_prob < 0.0 or neg_keep_prob > 1.0:
        raise ValueError("--neg-keep-prob must be between 0 and 1")

    chrom_offset = sum(ord(c) for c in str(chrom))
    rng = random.Random(seed + chrom_offset)

    kept = []
    for row in rows:
        if row["LABEL"] == 1 or rng.random() < neg_keep_prob:
            kept.append(row)

    return kept


def make_position_bed(rows, chrom, tmpdir, out_tsv_gz):
    positions = sorted({row["POS"] for row in rows})

    if tmpdir is None:
        bed_path = Path(str(out_tsv_gz) + ".pos.bed")
    else:
        bed_path = Path(tmpdir) / (Path(out_tsv_gz).name + ".pos.bed")

    with bed_path.open("w") as bed:
        for pos1 in positions:
            bed.write(f"{chrom}\t{pos1 - 1}\t{pos1}\n")

    return bed_path


def load_reference_sequence(ref_path, chrom):
    fasta = pysam.FastaFile(ref_path)
    seq = fasta.fetch(chrom).upper()
    fasta.close()
    return seq


def ctx5(seq, pos1):
    i = pos1 - 1
    start = max(0, i - 2)
    end = min(len(seq), i + 3)

    fragment = seq[start:end]

    if len(fragment) < 5:
        if start == 0:
            fragment = ("N" * (5 - len(fragment))) + fragment
        else:
            fragment = fragment + ("N" * (5 - len(fragment)))

    return fragment


def gc11(seq, pos1):
    i = pos1 - 1
    start = max(0, i - 5)
    end = min(len(seq), i + 6)

    fragment = seq[start:end]

    if not fragment:
        return 0.0

    gc = sum(1 for base in fragment if base in ("G", "C"))
    return gc / len(fragment)


def homopolymer_len(seq, pos1):
    i = pos1 - 1

    if i < 0 or i >= len(seq):
        return 0

    base = seq[i]

    if base == "N":
        return 1

    left = i
    while left - 1 >= 0 and seq[left - 1] == base:
        left -= 1

    right = i
    while right + 1 < len(seq) and seq[right + 1] == base:
        right += 1

    return right - left + 1


def parse_bases_quals(bases, quals):
    q_index = 0

    ref_count = 0

    snp_counts = defaultdict(int)
    snp_fwd = defaultdict(int)
    snp_rev = defaultdict(int)
    snp_bq = defaultdict(list)

    indel_counts = defaultdict(int)
    indel_fwd = defaultdict(int)
    indel_rev = defaultdict(int)

    i = 0
    length = len(bases)

    while i < length:
        char = bases[i]

        # Start of read marker followed by mapping quality character
        if char == "^" and i + 1 < length:
            i += 2
            continue

        # End of read marker
        if char == "$":
            i += 1
            continue

        # Placeholder for deletion in subsequent pileup columns
        if char == "*":
            if q_index < len(quals):
                q_index += 1
            i += 1
            continue

        # Insertion or deletion annotation
        if char in "+-":
            match = INDEL_RE.match(bases, i)

            if not match:
                i += 1
                continue

            sign = match.group(1)
            indel_len = int(match.group(2))
            seq_start = match.end()
            seq = bases[seq_start : seq_start + indel_len]

            if len(seq) == indel_len:
                key = (sign, seq.upper())
                indel_counts[key] += 1

                if any(base.islower() for base in seq):
                    indel_rev[key] += 1
                else:
                    indel_fwd[key] += 1

            i = seq_start + indel_len
            continue

        base_quality = None
        if q_index < len(quals):
            base_quality = ord(quals[q_index]) - 33

        # Reference match
        if char in ".,":
            ref_count += 1

            if base_quality is not None:
                q_index += 1

            i += 1
            continue

        # SNP observation
        if char.upper() in ("A", "C", "G", "T", "N"):
            base = char.upper()

            snp_counts[base] += 1

            if char.islower():
                snp_rev[base] += 1
            else:
                snp_fwd[base] += 1

            if base_quality is not None:
                snp_bq[base].append(base_quality)
                q_index += 1

            i += 1
            continue

        if base_quality is not None:
            q_index += 1

        i += 1

    dp_observed = len(quals)

    return {
        "dp_observed": dp_observed,
        "ref_count": ref_count,
        "snp_counts": snp_counts,
        "snp_fwd": snp_fwd,
        "snp_rev": snp_rev,
        "snp_bq": snp_bq,
        "indel_counts": indel_counts,
        "indel_fwd": indel_fwd,
        "indel_rev": indel_rev,
    }


def run_mpileup(args, bed_path):
    mpileup_cmd = [
        "samtools",
        "mpileup",
        "-r",
        args.chrom,
        "-l",
        str(bed_path),
        "-f",
        args.ref,
        "-q",
        str(args.min_mapq),
        "-Q",
        str(args.min_bq),
        "-d",
        str(args.max_depth),
        args.bam,
    ]

    process = subprocess.Popen(
        mpileup_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    assert process.stdout is not None

    stats_by_pos = {}

    for line in process.stdout:
        parts = line.rstrip().split()

        if len(parts) < 6:
            continue

        pos1 = int(parts[1])
        ref_base = parts[2].upper()
        dp = int(parts[3])
        bases = parts[4]
        quals = parts[5]

        parsed = parse_bases_quals(bases, quals)

        stats_by_pos[pos1] = {
            "ref_base": ref_base,
            "dp": dp,
            **parsed,
        }

    process.stdout.close()

    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()

    if return_code != 0:
        sys.stderr.write(stderr)
        raise SystemExit(f"samtools mpileup failed with code {return_code}")

    return stats_by_pos


def extract_variant_features(row, stats_by_pos):
    pos1 = row["POS"]
    ref = row["REF"]
    alt = row["ALT"]

    stats = stats_by_pos.get(pos1)

    if stats is None:
        return {
            "DP": 0,
            "ALT_COUNT": 0,
            "REF_COUNT": 0,
            "VAF": 0.0,
            "ALT_FWD": 0,
            "ALT_REV": 0,
            "SB": 0.0,
            "ALT_BQ_MEAN": "",
        }

    dp = int(stats["dp"])
    ref_count = int(stats["ref_count"])

    alt_count = 0
    alt_fwd = 0
    alt_rev = 0
    alt_bq_mean = ""

    # SNP
    if len(ref) == 1 and len(alt) == 1:
        base = alt.upper()

        alt_count = int(stats["snp_counts"].get(base, 0))
        alt_fwd = int(stats["snp_fwd"].get(base, 0))
        alt_rev = int(stats["snp_rev"].get(base, 0))

        base_qualities = stats["snp_bq"].get(base, [])
        if base_qualities:
            alt_bq_mean = sum(base_qualities) / len(base_qualities)

    # Insertion
    elif len(ref) == 1 and len(alt) > 1:
        inserted = alt[1:].upper()
        key = ("+", inserted)

        alt_count = int(stats["indel_counts"].get(key, 0))
        alt_fwd = int(stats["indel_fwd"].get(key, 0))
        alt_rev = int(stats["indel_rev"].get(key, 0))

    # Deletion
    elif len(ref) > 1 and len(alt) == 1:
        deleted = ref[1:].upper()
        key = ("-", deleted)

        alt_count = int(stats["indel_counts"].get(key, 0))
        alt_fwd = int(stats["indel_fwd"].get(key, 0))
        alt_rev = int(stats["indel_rev"].get(key, 0))

    vaf = alt_count / dp if dp > 0 else 0.0

    strand_total = alt_fwd + alt_rev
    sb = abs(alt_fwd - alt_rev) / strand_total if strand_total > 0 else 0.0

    return {
        "DP": dp,
        "ALT_COUNT": alt_count,
        "REF_COUNT": ref_count,
        "VAF": vaf,
        "ALT_FWD": alt_fwd,
        "ALT_REV": alt_rev,
        "SB": sb,
        "ALT_BQ_MEAN": alt_bq_mean,
    }


def write_features(args, rows, seq, stats_by_pos):
    out_path = Path(args.out_tsv_gz)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    tmp_out = Path(str(out_path) + ".tmp")

    with gzip.open(tmp_out, "wt") as out:
        out.write("\t".join(HEADER) + "\n")

        for row in rows:
            features = extract_variant_features(row, stats_by_pos)

            alt_bq_mean = features["ALT_BQ_MEAN"]

            out.write(
                f"{row['CHROM']}\t"
                f"{row['POS']}\t"
                f"{row['REF']}\t"
                f"{row['ALT']}\t"
                f"{row['TYPE']}\t"
                f"{row['SUPPORT']}\t"
                f"{row['LABEL']}\t"
                f"{features['DP']}\t"
                f"{features['ALT_COUNT']}\t"
                f"{features['REF_COUNT']}\t"
                f"{features['VAF']:.6f}\t"
                f"{features['ALT_FWD']}\t"
                f"{features['ALT_REV']}\t"
                f"{features['SB']:.6f}\t"
                f"{alt_bq_mean}\t"
                f"{ctx5(seq, row['POS'])}\t"
                f"{homopolymer_len(seq, row['POS'])}\t"
                f"{gc11(seq, row['POS']):.6f}\n"
            )

    os.replace(tmp_out, out_path)


def main():
    args = parse_args()

    if not Path(args.bam).is_file():
        raise SystemExit(f"ERROR: missing BAM: {args.bam}")

    if not Path(args.ref).is_file():
        raise SystemExit(f"ERROR: missing reference FASTA: {args.ref}")

    if not Path(args.labeled_tsv_gz).is_file():
        raise SystemExit(f"ERROR: missing labeled TSV: {args.labeled_tsv_gz}")

    rows = read_labeled_rows(args)
    rows = downsample_negatives(rows, args.neg_keep_prob, args.seed, args.chrom)

    if not rows:
        write_header_only(args.out_tsv_gz)
        print(f"No rows for {args.chrom}. Wrote header only.", file=sys.stderr)
        return

    bed_path = make_position_bed(rows, args.chrom, args.tmpdir, args.out_tsv_gz)

    try:
        seq = load_reference_sequence(args.ref, args.chrom)
        stats_by_pos = run_mpileup(args, bed_path)
        write_features(args, rows, seq, stats_by_pos)

    finally:
        bed_path.unlink(missing_ok=True)

    print(f"Wrote {args.out_tsv_gz}", file=sys.stderr)


if __name__ == "__main__":
    main()
