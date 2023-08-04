#!/usr/bin/env python3
# Gets statistics for read events.

import pdb

import pandas

# FIXME serialize this?
print('[reading in events]')
events = pandas.read_csv('read_event_table_5000.csv.gz')   # , nrows=100)
print('[computing stats]')
event_stats = events.describe()
print('[writing stats]')
event_stats.to_csv('read_event_stats.csv')

