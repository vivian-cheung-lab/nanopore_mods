#!/usr/bin/env python3
# Computes PCA of current and dwell time.

import pdb

import numpy as np
import pandas
import sklearn.decomposition
import sklearn.preprocessing

# read in table
event_table = pandas.read_csv('read_event_table.csv', index_col=0)
# convert to numpy
X = event_table.iloc[:,1:].to_numpy()

def standardize(x):
    """Standardizes a data matrix."""
    scaler = sklearn.preprocessing.StandardScaler()
    scaler.fit(x)
    return scaler.transform(x)

X = standardize(X)
pca = sklearn.decomposition.PCA(n_components=3)
Y = pca.fit_transform(X)


pdb.set_trace()

# ??? or standardize rows?

