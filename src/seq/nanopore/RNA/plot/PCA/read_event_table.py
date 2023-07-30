#!/usr/bin/env python3
# Writes out a table of event current and dwell time (e.g. for PCA).

import pdb
import sys

import numpy as np
import pandas



snakemake_dir = '../../polish/f5c/IVT/mods/Snakemake'
matrix_dir = f'{snakemake_dir}/event_matrix/chr19:45406985-45408892/'

sample_table = pandas.read_csv(f'{snakemake_dir}/IVT_samples.csv')
sample_table['name_and_concentration'] = (sample_table.short_name
    + ' ' + sample_table.mod_concentration)

def read_event_matrices():
    """Reads in the events."""
    print('[reading events]')
    matrices = dict()
    for r in sample_table[:10].itertuples():
        print('.', end='')
        event_data = np.load(f'{matrix_dir}/{r.FASTQ_file_base_name}.npz')
        m = event_data['event_stats']
        matrices[r.name_and_concentration] = m
        sys.stdout.flush()
    print()



def get_count_of_nonmissing_ranges(event_stats):
    pass    


m = read_event_matrices()


