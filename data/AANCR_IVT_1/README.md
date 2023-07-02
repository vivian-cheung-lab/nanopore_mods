# Sample data for nanopore_mod pipeline

This directory contains some sample data for running this pipeline.

- `AANCR_IVT.fa`: the expected sequence of the RNA. In this case,
this was the sequence of the clone used for IVT. (Alternatively,
this could be, e.g., the 
[human GRCh38 reference sequence](
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz)).

- `regions.bed`: the regions to include

- `mop_preprocess/`: this directory should contain the output of
the [Master of Pores (version 2) pipeline](
https://biocorecrg.github.io/MOP2/docs/about.html).
Because this data is large, it's distributed separately as a
`.zip` archive (which will need to be unpacked).
