# nanopore_mods

Code associated with the paper

*Nanopore-based direct sequencing of RNA transcripts with ten different modified nucleotides reveals gaps in existing technology.*
Joshua Burdick, Annelise Comai, Alan Bruzel, Vivian Cheung

### Master of Pores 2

We used the preprocessing and modification calling
workflows from
[Master of Pores 2](https://github.com/biocorecrg/MOP2).
For running those steps, please see
[usage of MOP2](MOP2).


### Snakemake

Several of the steps were automated using
[Snakemake](https://snakemake.readthedocs.io/en/stable/).
To run these, use something like:

```
cd workflow
snakemake --cores 8
```

