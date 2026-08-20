#!/usr/bin/env python3
# Writes a reference sequence for AANCR.
# This is basically just chromosome 19, with a few bases changed
# (based on the sequencing of the clone).
#
# The overall bounds of the clone are chr19:45406985-45408892 .

import os
import subprocess

import pyfaidx

output_fasta = 'chr19_AANCR_clone_only_230411.fa'

# download chromosome 19 from UCSC
if True:
    subprocess.run(['wget',
        'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz'])
    subprocess.run(['gunzip', 'chr19.fa.gz'])
    subprocess.run(['cp', 'chr19.fa', output_fasta])

fa = pyfaidx.Fasta(output_fasta, mutable=True)

print(fa['chr19'][45407786:45407790])
# change 45407788 (in 1-based coordinates) G>A
# !!! note that this alters the file on-disk!
fa['chr19'][45407787] = 'a'
print(fa['chr19'][45407786:45407790])
print()
# also change this base
print(fa['chr19'][45408833:45408837])
fa['chr19'][45408835] = 'G'
print(fa['chr19'][45408833:45408837])

# also, mask the rest of chr19 (outside of AANCR) to N's
def mask_region(start, end):
    fa['chr19'][start:end] = 'N' * (end-start)
# considering AANCR bounds to be chr19:45406985-45408892
mask_region(0, 45406984)
mask_region(45408893, len(fa['chr19']))

