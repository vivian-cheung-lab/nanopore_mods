#!/usr/bin/env python3
# Gets statistics for read events.

import pdb

import pandas
import scipy.stats

# FIXME serialize this?
print('[reading in events]')
# for practice, using 500 rows of previous table; this should include > 1 sample
# events = pandas.read_csv('read_event_table.csv', nrows=500)
events = pandas.read_csv('read_event_table_bases_5000reads.csv.gz')
print('[computing stats]')
event_stats = events.describe()
print('[writing stats]')
# event_stats.to_csv('read_event_stats.csv')

events_by_sample = events.set_index('sample name').iloc[:,1:]
sample_names = events_by_sample.index.unique()

def column_ANOVA(column_name):
    """Does a one-way ANOVA of a column."""
    print(column_name, end='                                  \r')
    # first, reshape numbers into a list of lists
    events_1 = events_by_sample[column_name]
    x = [events_1[s].tolist() for s in sample_names]
    # then, pass in the list of lists as args
    return scipy.stats.f_oneway(*x)

column_ANOVAs = [column_ANOVA(c) for c in events_by_sample.columns]
column_ANOVAs = pandas.DataFrame(column_ANOVAs,
        index=events_by_sample.columns)
column_ANOVAs.to_csv('read_event_stats_ANOVA_bases_5000reads.csv.gz')

