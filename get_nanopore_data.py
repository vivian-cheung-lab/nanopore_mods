#!/usr/bin/env python3
# Downloads FAST5 files, and gets FASTQ files from them.

import glob
import gzip
import os
import pdb
import subprocess

import pandas

from ont_fast5_api.fast5_interface import get_fast5_file
from pysradb import SRAweb

# Where to write output.
output_dir = "IVT_AANCR/"

os.makedirs(output_dir, exist_ok=True)

# write out just the clone sequence on chr19
subprocess.run(["./chr19_AANCR_clone_only.py"])

# Get sample info from SRA.
sra_client = SRAweb()
sample_info = sra_client.sra_metadata("SRP439844", detailed=True)
sample_info["run_total_bases"] = sample_info["run_total_bases"].astype(int)
sample_info.to_csv(f"{output_dir}/SRP439844_sample_info.csv")
# Original sample names.
original_samples = pandas.read_csv("IVT_samples.csv")
# Tweak names to match those in the SRA annotation.
original_samples.short_name.replace({
        "Y": "pseudouridine"
    },
    inplace=True)
original_samples['mod_concentration'] = (
    original_samples['mod_concentration'].str.replace(":", " to "))
original_samples.mod_concentration.replace({
    "100%": "1",
    },
    inplace=True)
original_samples.columns = [
    "FASTQ_file_base_name",
    "modification_name",
    "modification_concentration"
]
sample_info = sample_info.merge(original_samples)

# for testing: start with smallest samples
sample_info = sample_info.sort_values("run_total_bases").head(2)

def format_file_name(s):
    """Makes a filename.

    Rewrites a sample name for a filename, avoiding some characters
    (such as ':' and '/') which might cause problems in filenames.
    """
    s = s.replace(':', '_to_')
    s = s.replace('w/', 'with')
    s = s.replace('%', 'pct')
    s = s.replace('(', '')
    s = s.replace(')', '')
    s = s.replace(' ', '_')
    return s

def write_fastq(fast5_filepath, f):
    """Writes reads from one FAST5 file to a file.

    fast5_filepath: name of the (basecalled) FAST5 file
    f: file-like object to write reads to
    Side effects: writes reads to f
    """
    # We only write "pass" reads, which have a mean qscore >= this.
    min_pass_qscore = 7
    with get_fast5_file(fast5_filepath, mode="r") as f5:
        for read in f5.get_reads():
            # get sequence and quality, as a FASTQ-format string
            fastq = read.handle['Analyses/Basecall_1D_000/BaseCalled_template']['Fastq'][()]
            fastq_str = fastq.tobytes().decode()
            f.write(fastq_str)

def write_all_fastq(fast5_dir, output_file):
    """Writes a FASTQ file of reads from all the FAST5 files.

    fast5_dir: directory containing FAST5 files
    output_file: .fastq.gz file to write
    Side effects: writes `output_file`
    """
    with gzip.open(output_file, "wt") as fastq_file:
        for fast5_file in glob.glob(f"{fast5_dir}/*.fast5"):
            print(f"writing FASTQ from {fast5_file}")
            write_fastq(fast5_file, fastq_file)

def write_sample(srr_id, library_name, output_name):
    """Writes files for one sample.
    
    srr_id: the SRR run ID
    library_name: name of the sample
    output_name: name of directory to write
    Side effects:
    - Downloads the FAST5 files.
    - Extracts the FASTQ reads from these.
    """
    # XXX The directory names aren't exactly what the Snakemake
    # workflow expects, but should be close enough.
    output_dir_1 = f"{output_dir}/{output_name}/{library_name}/reads/"
    os.makedirs(f"{output_dir_1}/fast5/", exist_ok=True)

    # Write FAST5 files. The URL for the data is e.g.:
    # https://sra-pub-src-1.s3.amazonaws.com/SRR25667111/hm5C_1_hm5C_to_2_C.tar.gz.1
    tar_gz_filename = format_file_name(library_name) + ".tar.gz"
    fast5_url = f"https://sra-pub-src-1.s3.amazonaws.com/{srr_id}/{tar_gz_filename}.1"
    print(f"fast5_url = {fast5_url}")
    subprocess.run(f"curl {fast5_url} | tar xvfz - --strip-components=1",
        shell=True, check=True, cwd=f"{output_dir_1}/fast5/")

    # Write FASTQ files.    
    os.makedirs(f"{output_dir_1}/guppy_output/", exist_ok=True)    
    write_all_fastq(
        f"{output_dir_1}/fast5/",
        f"{output_dir_1}/guppy_output/pass.fastq.gz")

for _, row in sample_info.iterrows():
    write_sample(
        row["run_accession"],
        row["library_name"],
        row["FASTQ_file_base_name"]
    )
