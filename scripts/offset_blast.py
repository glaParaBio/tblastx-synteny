#!/usr/bin/env python3


import argparse
import subprocess as sp
import sys
import re


def header_dict(lst):
    if len(lst) != len(set(lst)):
        raise Exception("Duplicate column names")
    for x in ["qaccver", "qstart", "qend"]:
        if x not in lst:
            raise Exception(f'"{x}" not found in header')
    header = {}
    for i in range(0, len(lst)):
        header[lst[i]] = i
    return header


def parse_region(reg):
    parts = reg.split(":")
    if len(parts) < 2:
        raise Exception(
            f"Region format expected to be <chrom>:<start>-<end>. Got: {reg}"
        )
    ctg = ":".join(parts[0 : len(parts) - 1])
    start, end = parts[-1].split("-")
    return ctg, int(start), int(end)


parser = argparse.ArgumentParser(
    description="Offset the coordinates in blast output from a sequence chunk",
    formatter_class=lambda prog: argparse.RawDescriptionHelpFormatter(prog, width=80),
)

parser.add_argument(
    "blast",
    help="Blast output to process [%(default)s]",
    nargs="?",
    default="-",
)
parser.add_argument("--version", action="version", version="%(prog)s 0.1.0")

args = parser.parse_args()

if args.blast == "-":
    fin = sys.stdin
else:
    fin = open(args.blast)

header = {}

for line in fin:
    line = line.strip().split("\t")
    if not header:
        header = header_dict(line)
        print("\t".join(line))
        continue
    ctg, start, end = parse_region(line[header["qaccver"]])
    offset = start - 1
    assert offset >= 0
    line[header["qaccver"]] = ctg
    line[header["qstart"]] = str(int(line[header["qstart"]]) + offset)
    line[header["qend"]] = str(int(line[header["qend"]]) + offset)
    print("\t".join(line))

if fin != sys.stdin:
    fin.close()
