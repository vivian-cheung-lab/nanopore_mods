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

# range of 0-based positions to include (half-open)
# range_to_include = (10, 110)
range_to_include = (10, 210)

def read_event_matrices():
    """Reads in the events."""
    print('[reading events]')
    matrices = dict()
    for r in sample_table.itertuples():
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
    # get event matrices for unmodified samples
    unmodified = [matrices[s]
            for s in matrices.keys()
            if s.endswith(' none')]
    # stack the matrices
    x = np.concatenate(unmodified, axis=0)
    # get statistics
    r = {
            'n': np.sum(~np.isnan(x), axis=0),
            'mean': np.nanmean(x, axis=0),
            'std': np.nanstd(x, axis=0)
    }
    return r


matrices = read_event_matrices()
print('[getting unmodified stats]')
unmodified_stats = get_means_for_missing_values()
# column names
column_names = (
        [f'current {i+1}' for i in range(range_to_include[0], range_to_include[1])]
        + [f'time {i+1}' for i in range(range_to_include[0], range_to_include[1])])
# unmodified stats for those locations
unmodified_1 = (unmodified_stats['mean'][:,range_to_include[0]:range_to_include[1]].flatten())


def format_array(x):
    """Formats one array of event stats with those numbers."""
    # extract just the selected columns
    x1 = x[:,:,range_to_include[0]:range_to_include[1]]
    # sort by number of "present" numbers, decreasing
    num_present = np.sum(~np.isnan(x1[:,0,:]), axis=1)
    i = np.argsort(num_present)[-100:]
    # select the rows with the highest numbers
    x1 = x1[i,:,:]
    # reshape to put current first, then time
    # XXX there's probably a better way to do this
    x2 = np.concatenate([x1[:,0,:], x1[:,1,:]], axis=1)
    # fill in missing values with mean
    m1 = np.broadcast_to(unmodified_1, x2.shape)
    i = np.isnan(x2)
    x2[i] = m1[i]
    # round the numbers slightly
    x2 = np.round(x2, 6)
    return x2

def format_as_table(sample_name):
    """Formats one sample as a Pandas dataframe."""
    x = matrices[sample_name]
    x = format_array(x)
    x = pandas.DataFrame(x, columns=column_names)
    # XXX prepend sample name
    x['sample name'] = sample_name
    n = x.shape[1]
    j = [n-1] + list(range(n-1))
    x = x.iloc[:,j]
    return x

print('[formatting tables]')
event_tables = [format_as_table(s) for s in matrices.keys()]

print('[writing table]')
event_table = pandas.concat(event_tables)
event_table.to_csv('read_event_table.csv')

