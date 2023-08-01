#!/usr/bin/env python3
# Computes PCA of current and dwell time.

import pdb

import numpy as np
import pandas
import sklearn.decomposition
import sklearn.preprocessing

# read in table
event_table = pandas.read_csv('read_event_table.csv', index_col=0)
event_table.reset_index(drop=True, inplace=True)

# convert to numpy
X = event_table.iloc[:,1:].to_numpy()

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

# write out with columns standardized (as before)
write_PCA(standardize(X), 'PCA_columns_standardized.csv')
# ... and with rows standardized (treating current and time separately)
n = X.shape[1] // 2
X1 = np.concatenate([
    standardize(X[:,:n].T).T,
    standardize(X[:,n:(2*n)].T).T,
    ], axis=1)
write_PCA(X1, 'PCA_rows_standardized.csv')

