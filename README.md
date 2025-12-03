<!-- vim-markdown-toc GFM -->

* [Setup](#setup)
* [Run](#run)
* [Notes & Tips](#notes--tips)
* [Developer](#developer)

<!-- vim-markdown-toc -->

`tblastx-synteny` generates synteny files in
[paf](https://lh3.github.io/minimap2/minimap2.html#10) format using tblastx or
blastn or megablast. Compared to e.g. minimap2, `tblastx-synteny` is more sensitive in detecting regions
of synteny bewteen distant species since it uses blast (at the cost of computational time and scalability).

To visualize paf alignments you can use [Apollo3](https://apollo.jbrowse.org/)
(for Apollo see also this
[tutorial](https://github.com/glaParaBio/apollo-synteny-example?tab=readme-ov-file))

# Setup

Get the source code of this project using git or dowload the zip file from GitHub:

```
git clone https://github.com/glaParaBio/tblastx-synteny
cd tblastx-synteny
```

Install dependencies. We install with conda in a dedicated environment.
Alternatively, see dependencies in [requirements.txt](requirements.txt) and use
your favourite installation method.

```
conda create -n tblastx-synteny
conda activate tblastx-synteny
conda install --yes -n tblastx-synteny --file requirements.txt
```

Install:

```
pip install .
```

# Run

This may take a couple of minutes:

```
conda activate tblastx-synteny

tblastx-synteny \
    --subject test/data/PlasmoDB-68_PbergheiANKA_Genome.fasta \
    --query test/data/Pfalciparum3D7.small.fasta \
    --query-regions test/data/query_regions.bed \
    --directory out \
    --chunk-size 20000
```

Output will be in the argument passed to `-d/--directory`.

# Notes & Tips

* Blast, especially tblastx, is slow. Use `tblastx-synteny` on genomes larger
  than a few  tens of Mb is probably impractical. Consider using
  `-r/--query-regions` to restrict to regions of interest

* Use masked references to speed up blast

* tblastx may produce many short hits even when aligning similar sequences or
  even when aligning a sequence to itself. This may be due to blast masking low
  complexity amino-acid sequences that will not be aligned. To disable masking
  add to `--blast-args` (i.e. to tblastx) the option `-seg no`

* For tblastx the tags `pi` (percentage identical matches) and `ni` (number of
  identical matches) refer to the protein to protein alignment.

# Developer

Run tests:

```
./test/test.py
```
