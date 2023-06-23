# nanopore_mods

Code associated with the paper

*Nanopore-based direct sequencing of RNA transcripts with ten different modified nucleotides reveals gaps in existing technology.*
Joshua Burdick, Annelise Comai, Alan Bruzel, Vivian Cheung

## Steps to reproduce analysis

The following describes steps needed to reproduce the analysis.

### Initial basecalling

Bases were initially called using Guppy, on a GPU-equipped cluster,
using the script in `workflow/scrips/run_guppy_batch.sh`.

### Snakemake

Several of the steps were automated using
[Snakemake](https://snakemake.readthedocs.io/en/stable/).
To run these, use something like:

```
cd workflow
snakemake --cores 8
```

### Master of Pores 2

We also used the preprocessing and modification calling
workflows from
[Master of Pores 2](https://github.com/biocorecrg/MOP2).
For running those steps, please see
[usage of MOP2](MOP2).
