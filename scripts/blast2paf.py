#!/usr/bin/env python3

import argparse
import subprocess as sp
import sys
import re
import pandas


def prepare_fai_dict(fn):
    fai = {}
    with open(fn) as fin:
        for line in fin:
            line = line.strip().split("\t")
            assert line[0] not in fai
            fai[line[0]] = line[1]
    return fai


parser = argparse.ArgumentParser(
    description="Convert blast to paf format",
    formatter_class=lambda prog: argparse.RawDescriptionHelpFormatter(prog, width=80),
)

parser.add_argument(
    "--blast",
    "-b",
    help="Blast output to be reformatted",
    required=True,
)
parser.add_argument(
    "--query-fai",
    "-q",
    help="Fasta index of query",
    required=True,
)
parser.add_argument(
    "--blast-task",
    "-t",
    choices=["tblastx", "blastn", "megablast"],
    help="Blast task that produced the input",
    required=True,
)
parser.add_argument(
    "--subject-fai",
    "-s",
    help="Fasta index of subject",
    required=True,
)
parser.add_argument(
    "--min-pident",
    "-p",
    type=float,
    help="Keep alignments with at least this percent identity [%(default)s]",
    default=0,
)
parser.add_argument(
    "--max-evalue",
    "-e",
    type=float,
    help="Keep alignments with at most this e-value [%(default)s]",
    default=10,
)
parser.add_argument(
    "--min-length",
    "-l",
    type=int,
    help="Keep alignments with at least this nucleotide length [%(default)s]",
    default=0,
)

parser.add_argument("--header", "-H", action="store_true", help="Output header line")
parser.add_argument("--version", action="version", version="%(prog)s 0.1.0")

args = parser.parse_args()

query_fai = prepare_fai_dict(args.query_fai)
subject_fai = prepare_fai_dict(args.subject_fai)

out = pandas.read_csv(args.blast, sep="\t")
if args.blast_task == "tblastx":
    out["length"] = out["length"] * 3
out = out[
    (out.pident >= args.min_pident)
    & (out.evalue <= args.max_evalue)
    & (out.length >= args.min_length)
]
out.rename(columns={"#qaccver": "qaccver"}, inplace=True)
out = out.drop(["sstrand"], axis=1)

qstart = []
qend = []
sstart = []
send = []
qaccver_length = []
saccver_length = []
strand = []
for row in out.itertuples():
    qaccver_length.append(query_fai[row.qaccver])
    saccver_length.append(subject_fai[row.saccver])
    if ((row.qstart < row.qend) and (row.sstart < row.send)) or (
        (row.qstart > row.qend) and (row.sstart > row.send)
    ):
        strand.append("+")
    else:
        strand.append("-")

    if row.qstart < row.qend:
        qstart.append(row.qstart)
        qend.append(row.qend)
    else:
        qstart.append(row.qend)
        qend.append(row.qstart)
    if row.sstart < row.send:
        sstart.append(row.sstart)
        send.append(row.send)
    else:
        sstart.append(row.send)
        send.append(row.sstart)
out["qstart"] = [x - 1 for x in qstart]
out["sstart"] = [x - 1 for x in sstart]
out["qend"] = qend
out["send"] = send
out["qaccver_length"] = qaccver_length
out["saccver_length"] = saccver_length
out["strand"] = strand
out["mapq"] = 255
out["blast_pident"] = ["pi:f:%s" % x for x in out["pident"]]
# out["blast_nident"] = ["ni:i:%s" % x for x in out["nident"]]
out["blast_evalue"] = ["ev:f:%s" % x for x in out["evalue"]]
# out["nident"] = "*"

# Keep highest bitscore for alignments at the same coords
coords = [
    "qaccver",
    "saccver",
    "qstart",
    "qend",
    "sstart",
    "send",
    "strand",
]
best = out.groupby(coords, as_index=False).agg(bitscore=("bitscore", "max"))
out = out.merge(best)

out["blast_bitscore"] = ["bs:f:%s" % x for x in out["bitscore"]]
out = out.drop(["bitscore"], axis=1)
out = out.sort_values(["qaccver", "qstart", "qend", "saccver", "sstart", "send"])
out.rename(columns={"qaccver": "#qaccver"}, inplace=True)
out = out[
    [
        "#qaccver",
        "qaccver_length",
        "qstart",
        "qend",
        "strand",
        "saccver",
        "saccver_length",
        "sstart",
        "send",
        "nident",
        "length",
        "mapq",
        "blast_pident",
        "blast_nident",
        "blast_evalue",
        "blast_bitscore",
    ]
].drop_duplicates()
out.to_csv(sys.stdout, sep="\t", index=False, header=args.header)
