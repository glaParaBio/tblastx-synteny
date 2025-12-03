#!/usr/bin/env python3

import argparse
from pyfaidx import Fasta
import os
import sys

parser = argparse.ArgumentParser(
    description="Split fasta: one file per interval in bed file",
    formatter_class=lambda prog: argparse.RawDescriptionHelpFormatter(prog, width=80),
)

parser.add_argument(
    "fasta",
    help="Fasta file to be split",
)
parser.add_argument(
    "--bed",
    "-b",
    help="Bed file of intervals [%(default)s]",
    default="-",
)
parser.add_argument(
    "--outdir",
    "-o",
    help="Output directory [%(default)s]",
    default=".",
)

parser.add_argument(
    "--n-files-per-dir",
    "-n",
    type=int,
    help="Not implemented yet: If > 0, write files in sub-dirs containing at most these many files each [%(default)s]",
    default=0,
)

args = parser.parse_args()


if args.bed == "-":
    fin = sys.stdin
else:
    fin = open(args.bed)

fasta = Fasta(args.fasta)

for line in fin:
    if line.startswith("#"):
        continue
    chrom, start, end = line.split("\t")[0:3]
    start = int(start)
    end = int(end)
    pstart = "%09d" % (start + 1)
    pend = "%09d" % end
    name = f"{chrom}:{pstart}-{pend}"

    with open(os.path.join(args.outdir, f"{name}.fa"), "w") as fout:
        fout.write(f">{name}\n")
        fout.write(fasta[chrom][start:end].seq + "\n")
