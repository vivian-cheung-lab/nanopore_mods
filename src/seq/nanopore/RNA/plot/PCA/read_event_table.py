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
    for r in sample_table[:3].itertuples():
        print('.', end='')
        event_data = np.load(f'{matrix_dir}/{r.FASTQ_file_base_name}.npz')
        m = event_data['event_stats']
        matrices[r.name_and_concentration] = m
        sys.stdout.flush()
    print()
    return matrices


def get_count_of_nonmissing_ranges(x):
    """Gets counts of non-missing ranges.

    FIXME not currently working
    """
    # we just use whether current is present
    # (as time is present at same sites)
    current_present = ~ np.isnan(x[:,0,:])
    # get mean of cumulative sum, up to each b ase
    sum_of_bases_to_left = np.cumsum(current_present, axis=1)
    # compute means of columns of this
    mean_left = sum_of_bases_to_left.mean(axis=0)
    # compute the means over the spans
    n = mean_left.shape[0]
    span_mean = np.zeros((n,n))
    for i in range(n-1):
        for j in range(i+1, n):
            span_mean[i,j] = (mean_left[j] - mean_left[i]) / (j-i)
    return span_mean

def get_means_for_missing_values():
    """Gets means at each position to fill in missing values.
    
    This will be current and time at each position, from the
    unmodified samples.
    """
    unmodified_samples = [matrices[s]
            for s in matrices.keys()
            if s.endswith(' none')]
    pdb.set_trace()





matrices = read_event_matrices()

present_count = [get_count_of_nonmissing_ranges(m) for m in matrices.values()]


pdb.set_trace()



