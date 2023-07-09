#!/usr/bin/env python3

import os
import pdb

import numpy as np
import pandas
import scipy.stats

input_dir = 'event_stats/'

# read in files
stats = [pandas.read_csv(input_dir + f) for f in os.listdir(input_dir)]
stats = pandas.concat(stats)
# stats = stats[:6]   # practice

# add on columns for stats
if False:
    stats['t'] = None
    stats['p'] = None
    stats['Adjusted p'] = None

mods = stats['Modification'].drop_duplicates().tolist()

stats1 = stats.set_index(['Modification','Sample'])

def get_stats(row):
    (mod, sample) = row
    if sample=='Unmodified':
        return pandas.DataFrame([{'t': None, 'p': None, 'Adjusted p':None}])
    a = stats1.loc[mod].loc['Modified']
    b = stats1.loc[mod].loc['Unmodified']
    r = scipy.stats.ttest_ind_from_stats(
        a['mean'], a['std'], a['count'],
        b['mean'], b['std'], b['count'])
    return pandas.DataFrame([{
        't': r.statistic,
        'p': r.pvalue,
        'Adjusted p': len(mods) * r.pvalue}])

p_stats = stats[['Modification', 'Sample']].apply(get_stats, axis=1)
p_stats = pandas.concat(p_stats.tolist())
stats = stats.reset_index()
p_stats = p_stats.reset_index()
p_stats = pandas.concat([stats, p_stats], axis=1)

p_stats.to_excel('Dwell time stats.xlsx')

# stats = pandas.DataFrame(stats, columns=['t', 'p', 'Adjusted p'])

