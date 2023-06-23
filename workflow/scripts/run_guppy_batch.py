#!/usr/bin/env python3
# Runs Guppy on a series of files.
# Assumed to be run in a dir containing many
# dirs. of MinKNOW output.
# Intended to be run as a SLURM job, by ./slurm_run_guppy.sh .
# FIXME

import glob, os, subprocess

# the original set of directories
dirs1 = [os.path.abspath(d + '/../')
        for d in glob.glob('*/*/*/fast5/')]
# fast5 files sometimes seem to be up a directory
dirs2 = [os.path.abspath(d + '/../')
        for d in glob.glob('*/*/fast5/')]
# for now, include only one set of dirs
dirs = dirs1

for d in dirs:
    # skip cases which have run already
    if os.path.exists(d + '/guppy_output'):
        continue
    print('Running Guppy in dir:', flush=True)
    print(d, flush=True)
    subprocess.run(['guppy_basecaller', '-i', 'fast5',
        '-s', 'guppy_output', '-c', 'rna_r9.4.1_70bps_hac.cfg',
        '-x', 'cuda:all', '--moves_out', '--post_out', '--bam_out',
        '--compress_fastq', '--records_per_fastq', '0',
        '--align_ref', '/home/jtburd/seq/align/index/minimap2/hg19_male.mmi',
        '--align_type', 'full', '--num_alignment_threads', '16'],
        cwd=d)
    print(flush=True)
