#!/usr/bin/env python3
# Reshapes PC scores, and/or prints out statistics for them.

import pdb

import pandas

PC_scores = pandas.read_table('PCscores.txt')

PC_score_mean = PC_scores.groupby('Label')[['PC1','PC2','PC3']].mean()

PC_score_mean.round(4).to_csv('PC_scores_mean.csv')

