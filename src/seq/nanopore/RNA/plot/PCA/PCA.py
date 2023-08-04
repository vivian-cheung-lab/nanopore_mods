#!/usr/bin/env python3
# Computes PCA of current and dwell time.

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
    pca = sklearn.decomposition.PCA(n_components=3)
    Y = pca.fit_transform(X)
    Y = pandas.DataFrame(Y, columns = ['PC1', 'PC2', 'PC3'])
    Y = np.round(Y, 5)
    r = pandas.concat([event_table['sample name'], pandas.DataFrame(Y)], axis=1)
    r.to_csv(output_filename)

# table with entire region, up to 5,000 reads per sample
event_table = pandas.read_csv('read_event_table_5000.csv.gz')         # , nrows=20000)
# read in ANOVA stats, for filtering informative events
column_ANOVAs = pandas.read_csv('read_event_stats_ANOVA.csv.gz', index_col=0)
column_ANOVAs['measurement'] = [s.split()[0] for s in column_ANOVAs.index]
column_ANOVAs['base'] = [s.split()[1] for s in column_ANOVAs.index]
# get columns with f-score above the 75% cutoff of f-score for each of these
ANOVA_stat_stats = column_ANOVAs.groupby('measurement')['statistic'].describe()
cutoffs = ANOVA_stat_stats['75%']
columns_above_cutoff = pandas.concat([
    column_ANOVAs[(column_ANOVAs.measurement=='current') & (column_ANOVAs.statistic >= cutoffs['current'])],
    column_ANOVAs[(column_ANOVAs.measurement=='time') & (column_ANOVAs.statistic >= cutoffs['time'])] ])
filtered_event_table = event_table[ columns_above_cutoff.index.tolist() ]
# write out PCA for this
filtered_event_table.reset_index(drop=True, inplace=True)
X = filtered_event_table.to_numpy()
write_PCA(standardize(X), 'PCA_5000_reads.csv.gz')

if True:
    # 500nt region table
    # read in table
    event_table = pandas.read_csv('read_event_table.csv', index_col=0)
    event_table.reset_index(drop=True, inplace=True)
    # convert to numpy
    X = event_table.iloc[:,1:].to_numpy()
    # write out with columns standardized (as before)
    write_PCA(standardize(X), 'PCA_columns_standardized.csv')

