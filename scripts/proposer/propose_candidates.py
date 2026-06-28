#!/usr/bin/env python3

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path


INDEL_RE = re.compile(r"([+-])(\d+)")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate permissive SNP/INDEL candidate sites from samtools mpileup."
    )

    parser.add_argument("--chrom", required=True, help="Chromosome/contig to process")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--ref", required=True, help="Reference FASTA")
    parser.add_argument("--out-vcf", required=True, help="Output candidate VCF")
    parser.add_argument("--out-bed", required=True, help="Output candidate BED")

    parser.add_argument("--min-mapq", type=int, default=1)
    parser.add_argument("--min-bq", type=int, default=1)
    parser.add_argument("--min-snp-support", type=int, default=1)
    parser.add_argument("--min-indel-support", type=int, default=1)
    parser.add_argument("--max-indel-len", type=int, default=80)
    parser.add_argument("--max-depth", type=int, default=7000)

    return parser.parse_args()


def parse_support(bases, max_indel_len):
    snp = {}
    indel = {}

    i = 0
    while i < len(bases):
        c = bases[i]

        # Start of read marker followed by mapping quality character
        if c == "^" and i + 1 < len(bases):
            i += 2
            continue

        # End of read marker
        if c == "$":
            i += 1
            continue

        # Placeholder for deletion in subsequent pileup columns
        if c == "*":
            i += 1
            continue

        # Insertion or deletion
        if c in "+-":
            match = INDEL_RE.match(bases, i)
            if not match:
                i += 1
                continue

            sign = match.group(1)
            length = int(match.group(2))
            seq_start = match.end()
            seq = bases[seq_start : seq_start + length]

            if len(seq) == length and length <= max_indel_len:
                key = (sign, seq.upper())
                indel[key] = indel.get(key, 0) + 1

            i = seq_start + length
            continue

        # Reference match
        if c in ".,":
            i += 1
            continue

        # SNP observation
        base = c.upper()
        if base in "ACGTN":
            snp[base] = snp.get(base, 0) + 1

        i += 1

    return snp, indel


def main():
    args = parse_args()

    bam = Path(args.bam)
    ref = Path(args.ref)
    out_vcf = Path(args.out_vcf)
    out_bed = Path(args.out_bed)

    if not bam.is_file():
        sys.exit(f"ERROR: missing BAM: {bam}")

    if not ref.is_file():
        sys.exit(f"ERROR: missing reference FASTA: {ref}")

    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    out_bed.parent.mkdir(parents=True, exist_ok=True)

    tmp_vcf = Path(str(out_vcf) + ".tmp")
    tmp_bed = Path(str(out_bed) + ".tmp")

    mpileup_cmd = [
        "samtools",
        "mpileup",
        "-r",
        args.chrom,
        "-f",
        str(ref),
        "-q",
        str(args.min_mapq),
        "-Q",
        str(args.min_bq),
        "-d",
        str(args.max_depth),
        str(bam),
    ]

    wrote = 0

    with tmp_vcf.open("w") as vcf, tmp_bed.open("w") as bed:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write(
            '##INFO=<ID=TYPE,Number=1,Type=String,Description="Candidate variant type: SNP/INS/DEL">\n'
        )
        vcf.write(
            '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Alternative allele support count from mpileup">\n'
        )
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        process = subprocess.Popen(
            mpileup_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        assert process.stdout is not None

        for line in process.stdout:
            parts = line.rstrip().split()

            if len(parts) < 5:
                continue

            chrom = parts[0]
            pos = int(parts[1])
            ref_base = parts[2].upper()
            bases = parts[4]

            snp, indel = parse_support(bases, args.max_indel_len)

            # SNP candidates
            for base, count in snp.items():
                if (
                    base != ref_base
                    and base in "ACGT"
                    and count >= args.min_snp_support
                ):
                    vcf.write(
                        f"{chrom}\t{pos}\t.\t{ref_base}\t{base}\t.\tPASS\t"
                        f"TYPE=SNP;SUPPORT={count}\n"
                    )
                    bed.write(f"{chrom}\t{pos - 1}\t{pos}\n")
                    wrote += 1

            # INDEL candidates
            for (sign, seq), count in indel.items():
                if count < args.min_indel_support:
                    continue

                if sign == "+":
                    vcf.write(
                        f"{chrom}\t{pos}\t.\t{ref_base}\t{ref_base}{seq}\t.\tPASS\t"
                        f"TYPE=INS;SUPPORT={count}\n"
                    )
                    bed.write(f"{chrom}\t{pos - 1}\t{pos}\n")
                    wrote += 1

                else:
                    vcf.write(
                        f"{chrom}\t{pos}\t.\t{ref_base}{seq}\t{ref_base}\t.\tPASS\t"
                        f"TYPE=DEL;SUPPORT={count}\n"
                    )
                    bed.write(f"{chrom}\t{pos - 1}\t{pos}\n")
                    wrote += 1

        process.stdout.close()
        stderr = process.stderr.read() if process.stderr is not None else ""
        return_code = process.wait()

    if return_code != 0:
        tmp_vcf.unlink(missing_ok=True)
        tmp_bed.unlink(missing_ok=True)
        sys.stderr.write(stderr)
        sys.exit(return_code)

    os.replace(tmp_vcf, out_vcf)
    os.replace(tmp_bed, out_bed)

    print(f"CANDIDATES_WRITTEN={wrote}", file=sys.stderr)


if __name__ == "__main__":
    main()
