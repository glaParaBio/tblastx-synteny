<!-- vim-markdown-toc GFM -->

* [Setup](#setup)
* [Run](#run)
* [Notes & Tips](#notes--tips)

<!-- vim-markdown-toc -->

# Setup

Install dependencies:

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

This may take a a few minutes:

```
tblastx-synteny \
    --subject test/data/PlasmoDB-68_PbergheiANKA_Genome.fasta \
    --query test/data/Pfalciparum3D7.small.fasta \
    --directory out \
    --chunk-size 20000
```

# Notes & Tips

* Use masked references

* tblastx may produce many short hits even when aligning similar sequences or
  even when aligning a sequence to itself. This may be due to blast masking low
  complexity amino-acid sequences that will not be aligned. To disable masking
  add to `--blast-args` (i.e. to tblastx) the option `-seg no`

* For tblastx the tags `pi` (percentage identical matches) and `ni` (number of
  identical mataches) refer to the protein to protein alignment.
