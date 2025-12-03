#!/usr/bin/env python3

import argparse
import os
import subprocess as sp
import sys
from colorama import Fore, Style
from pathlib import Path


def main():
    SCRIPT_DIR = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description="Create a synteny file in PAF format using tblastx",
        formatter_class=lambda prog: argparse.RawDescriptionHelpFormatter(
            prog, width=80
        ),
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
        help="Output directory",
        required=True,
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
        default="tblastx",
        choices=["tblastx", "blastn", "megablast"],
    )
    _default = "-evalue 0.01 -max_target_seqs 10"
    parser.add_argument(
        "--blast-args",
        "-b",
        help=f"A single string of further arguments to blast [{_default}]",
        default=[_default],
        nargs=1,
    )
    _default = "--printshellcmds"
    parser.add_argument(
        "--smk-args",
        "-a",
        help=f"A single string of further arguments to snakemake. E.g. '--keep-going --notemp' [{_default}]",
        default=[_default],
        nargs=1,
    )

    parser.add_argument("--version", "-v", action="version", version="%(prog)s 0.1.0")

    args = parser.parse_args()

    query = os.path.abspath(args.query)
    subject = os.path.abspath(args.subject)
    dry_run = "--dry-run" if args.dry_run else ""
    smk_args = " ".join(args.smk_args)
    blast_args = " ".join(args.blast_args).strip()

    if blast_args:
        blast_args = f"blast_args='{blast_args}'"

    snakecmd = f"""snakemake \
        --jobs {args.jobs} \
        -s {SCRIPT_DIR}/workflows/blast.smk \
        --directory {args.directory} \
        {dry_run} \
        --config query={query} \
                 subject={subject} \
                 chunk_size={args.chunk_size} \
                 task={args.task} \
                 {blast_args} \
                 cwd={SCRIPT_DIR} \
        {smk_args}
    """

    sys.stderr.write(Fore.GREEN + snakecmd + Style.RESET_ALL)

    p = sp.Popen(
        snakecmd,
        shell=True,
        stdout=sp.PIPE,
        stderr=sp.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True,
    )

    for line in iter(p.stdout.readline, ""):
        sys.stderr.write(Fore.YELLOW + line.strip() + "\n" + Style.RESET_ALL)
    p.stdout.close()
    p.wait()

    if p.returncode != 0:
        sys.stderr.write(Fore.RED + f"Error executing {snakecmd}\n" + Style.RESET_ALL)
        sys.exit(p.returncode)


if __name__ == "__main__":
    if sys.argv[0].endswith(".exe"):
        sys.argv[0] = sys.argv[0][:-4]
    sys.exit(main())
