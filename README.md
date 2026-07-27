
# nanopore_mods

This is code associated with the paper

Burdick JT, Comai A, Bruzel A, Sun G, Dedon PC, Cheung VG. Nanopore-based direct sequencing of RNA transcripts with 10 different modified nucleotides reveals gaps in existing technology. G3 (Bethesda). 2023 Nov 1;13(11):jkad200. doi: 10.1093/g3journal/jkad200. PMID: (37655917)[https://pubmed.ncbi.nlm.nih.gov/37655917/]; PMCID: PMC10627276.

## Running nanopore analysis

This describes how to run the nanopore modification analysis, starting from Guppy's
base-called output. Commands to run at the shell prompt are shown in this font:

```
ls
```

1. You will need an x86 Linux system with
- at least 80 GB of disk space
- at least 32 GB RAM
- at least 8 processors are recommended
- a recent version of Conda installed

2. Install a version of conda. This may already be installed on your
system. If not, Miniconda is a good alternative. To install it,
go to the [Miniconda installer page](https://docs.conda.io/en/latest/miniconda.html),
and follow the instructions there.

3. Clone this repository, into a directory with enough disk space.

```
git clone https://github.com/vivian-cheung-lab/nanopore_mods.git
```

4. Create the Conda environment, and activate it.

```
cd nanopore_mods
conda env create -f nanopore_mods_1.yaml
conda activate nanopore_mods_1
```

5. Download and unpack the data. (This is 37 GB of data, and will take at least six minutes.
Half an hour may be more realistic, depending on your network connection.)

There are two alternative ways to do this.

- One way is to download from the cloud.  This command should be run in the
`nanopore_mods` directory (like the previous command).

```
curl https://storage.googleapis.com/nanopore_mods_test_data/IVT_AANCR_data_20230708.tar.gz | tar xvfz -
```

- Another alternative is to download the data from
the (SRA)[https://www.ncbi.nlm.nih.gov/sra].

To do this, from the `nanopore_mods` directory, run:
```

```

5. Run the 
[Snakemake](https://snakemake.readthedocs.io/en/stable/)
workflow.

The `--cores 8` argument indicates that 8 cores should be used.

```
cd src/seq/nanopore/RNA/polish/f5c/IVT/mods/Snakemake
snakemake --cores 8
```

This will write various output files in the directory tree.

### Master of Pores 2

In later versions of this code, we used the preprocessing and modification calling
workflows from
[Master of Pores 2](https://github.com/biocorecrg/MOP2).
For running those steps, please see
[usage of MOP2](MOP2).

