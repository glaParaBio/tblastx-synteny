import pandas
import glob

QUERY = config["query"]
SUBJECT = config["subject"]
BLAST_FMT = [
    "qaccver",
    "saccver",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "nident",
    "positive",
    "slen",
    "sstrand",
]

if config["task"] == "tblastx":
    task = "tblastx"
elif config["task"] == "blastn":
    task = "blastn -task blastn"
elif config["task"] == "megablast":
    task = "blastn -task megablast"
else:
    raise Exception(f"Invalid task for blast: {config['task']}")


qn = re.sub("\\.fa$|\\.fasta$|\\.fna", "", os.path.basename(QUERY), flags=re.I)
sn = re.sub("\\.fa$|\\.fasta$|\\.fna", "", os.path.basename(SUBJECT), flags=re.I)


localrules:
    all,


rule all:
    input:
        f"{qn}_vs_{sn}.paf",


rule makeblastdb:
    input:
        fa=SUBJECT,
    output:
        db=f"{SUBJECT}.nin",
    shell:
        r"""
        makeblastdb -in {input.fa} -dbtype nucl
        """


rule faidx:
    input:
        fa="{fasta}",
    output:
        fai="{fasta}.fai",
    shell:
        r"""
        samtools faidx {input.fa}
        """


checkpoint split_query:
    input:
        fa=QUERY,
        fai=f"{QUERY}.fai",
    output:
        outdir=directory("split_tmp"),
    params:
        chunk_size=config["chunk_size"],
        cwd=config["cwd"],
    shell:
        r"""
        mkdir {output.outdir}
        bedtools makewindows -g {input.fai} -w {params.chunk_size} \
        | {params.cwd}/scripts/split_fasta.py --bed - --outdir {output.outdir} {input.fa}
        """


rule blast:
    input:
        query="split_tmp/{window}.fa",
        db=f"{SUBJECT}.nin",
        subject=SUBJECT,
    output:
        tmp=temp("{blast}/{window}.tmp"),
        out="{blast}/{window}.out",
    params:
        fmt=" ".join(BLAST_FMT),
        task=task,
        blast_args="" if "blast_args" not in config else config["blast_args"],
        cwd=config["cwd"],
    shell:
        r"""
        echo "{params.fmt}" | tr ' ' '\t' > {output.tmp}
        {params.task} -query {input.query} \
           -db {input.subject} \
           {params.blast_args} \
           -outfmt "6 {params.fmt}" >> {output.tmp}
        {params.cwd}/scripts/offset_blast.py {output.tmp} > {output.out}
        """


rule blast_to_paf:
    input:
        out="%s/{window}.out" % config["task"],
        query_fai=f"{QUERY}.fai",
        subject_fai=f"{SUBJECT}.fai",
    output:
        paf=temp("paf/{window}.paf"),
    params:
        cwd=config["cwd"],
        task=config["task"],
    shell:
        r"""
        {params.cwd}/scripts/blast2paf.py \
                --blast {input.out} \
                --query-fai {input.query_fai} \
                --subject-fai {input.subject_fai} \
                --blast-task {params.task} \
                --header \
                --min-pident 0 \
                --min-length 0 \
                --max-evalue 1000 > {output.paf}
        """


def aggregate_paf(wc):
    checkpoint_output = checkpoints.split_query.get().output.outdir
    outfiles = [
        os.path.basename(x) for x in glob.glob(os.path.join(checkpoint_output, "*.fa"))
    ]
    windows = sorted([re.sub("\\.fa$", "", x) for x in outfiles])
    return expand("paf/{window}.paf", window=windows)


rule cat_paf:
    input:
        paf=aggregate_paf,
    output:
        paf=f"{qn}_vs_{sn}.paf",
    params:
        cwd=config["cwd"],
    shell:
        r"""
        head -n 1 {input.paf[0]} > {output.paf}
        cat {input.paf} | grep -v '^#' >> {output.paf}
        """
