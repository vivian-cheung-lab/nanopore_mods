#!/usr/bin/env python3
# Computes PCA of current and dwell time.

import os
import pdb

import numpy as np
import pandas
import sklearn.decomposition
import sklearn.preprocessing

def standardize(x):
    """Standardizes the columns of a data matrix."""
    scaler = sklearn.preprocessing.StandardScaler()
    scaler.fit(x)
    return scaler.transform(x)

def write_PCA(X, output_filename):
    """Computes and writes out a PCA."""
    # using the 'arpack' solver seems to be running faster than the
    # 'randomized' solver (which the default 'auto' setting picks)
    pca = sklearn.decomposition.PCA(n_components=3, svd_solver='arpack')
    Y = pca.fit_transform(X)
    Y = pandas.DataFrame(Y, columns = ['PC1', 'PC2', 'PC3'])
    Y = np.round(Y, 5)
    r = pandas.concat([event_table['sample name'], pandas.DataFrame(Y)], axis=1)
    r.to_csv(output_filename)

# table with entire region, up to 5,000 reads per sample
event_table = pandas.read_csv('read_event_table_bases_5000reads.csv.gz')     # , nrows=20000)
# read in ANOVA stats, for filtering informative events
column_ANOVAs = pandas.read_csv('read_event_stats_ANOVA_bases_5000reads.csv.gz', index_col=0)
column_ANOVAs['measurement'] = [s.split()[0] for s in column_ANOVAs.index]
column_ANOVAs['base'] = [s.split()[1] for s in column_ANOVAs.index]

def write_PCA_at_cutoff(quantile_cutoff):
    """Writes out a PCA using the most-variable numbers for some quantiles.

    quantile_cutoff: the fraction of most-variable

    Side effects: writes a CSV file of PC scores
    """
    print(f'[writing PCA with quantile cutoff of {quantile_cutoff}]')
    output_dir = 'PCA_bases_5000reads/'
    os.makedirs(output_dir, exist_ok=True)
    # find cutoff
    current_cutoff = np.quantile(
            column_ANOVAs[column_ANOVAs.measurement=='current'].statistic.to_numpy(),
            1 - quantile_cutoff)
    time_cutoff = np.quantile(
            column_ANOVAs[column_ANOVAs.measurement=='time'].statistic.to_numpy(),
            1 - quantile_cutoff)
    # extract columns with f-score above that cutoff
    columns_above_cutoff = pandas.concat([
        column_ANOVAs[(column_ANOVAs.measurement=='current') & (column_ANOVAs.statistic >= current_cutoff)],
        column_ANOVAs[(column_ANOVAs.measurement=='time') & (column_ANOVAs.statistic >= time_cutoff)] ])
    filtered_event_table = event_table[ columns_above_cutoff.index.tolist() ]
    # write out PCA for this
    filtered_event_table.reset_index(drop=True, inplace=True)
    X = filtered_event_table.to_numpy()
    write_PCA(standardize(X), f'{output_dir}/PCA_{quantile_cutoff}_quantile.csv.gz')

for qc in [.01, .05, .1, .25, .5]:
    write_PCA_at_cutoff(qc)

