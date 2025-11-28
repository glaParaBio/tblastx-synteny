#!/usr/bin/env python3

import argparse
import os
import subprocess as sp
import sys

parser = argparse.ArgumentParser(
    description="Create a synteny file in PAF format using tblastx",
    formatter_class=lambda prog: argparse.RawDescriptionHelpFormatter(prog, width=80),
)

parser.add_argument(
    "--query",
    "-q",
    help="Fasta file of query. This will be split in chunks",
    required=True,
)
parser.add_argument(
    "--subject",
    "-s",
    help="Fasta file of subject. This will be indexed by blast and index files will be in the same directory as the input file",
    required=True,
)
parser.add_argument(
    "--chunk-size",
    "-S",
    help="Split query in sequences of this size [%(default)s]",
    default=20000,
    type=int,
)
parser.add_argument(
    "--directory",
    "-d",
    help="Output directory [%(default)s]",
    default="synteny-output",
)
parser.add_argument(
    "--jobs",
    "-j",
    help="Number of jobs to run in parallel [%(default)s]",
    default=5,
    type=int,
)
parser.add_argument(
    "--dry-run",
    "-n",
    action="store_true",
    help="Show the commands that would be executed and exit",
)
parser.add_argument(
    "--task",
    "-t",
    help="Blast task to run [%(default)s]",
    default='tblastx',
    choices=['tblastx', 'blastn', 'megablast'],
)


parser.add_argument("--version", "-v", action="version", version="%(prog)s 0.1.0")

args = parser.parse_args()

query = os.path.abspath(args.query)
subject = os.path.abspath(args.subject)
dry_run = "--dry-run" if args.dry_run else ""
cwd = os.getcwd()

snakecmd = f"""
snakemake --printshellcmds \
    --jobs {args.jobs} \
    -s workflows/blast.smk \
    --directory {args.directory} \
    {dry_run} \
    --config query={query} \
             subject={subject} \
             chunk_size={args.chunk_size} \
             task={args.task} \
             cwd={cwd}
"""

print(snakecmd)

stderr = []
p = sp.Popen(snakecmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
while True:
    so = p.stdout.readline()
    se = p.stderr.readline()
    print(so, end="")
    print(se, end="")
    stderr.append(se.strip())
    if not so and not se:
        break
p.wait()

if p.returncode != 0:
    raise sp.CalledProcessError(
        returncode=p.returncode, cmd=snakecmd, output="\n".join(stderr)
    )
