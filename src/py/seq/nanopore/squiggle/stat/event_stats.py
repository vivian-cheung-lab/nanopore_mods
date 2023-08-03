# Event-level statistics, including:
# - getting statistics about events (e.g. mean, s.e.m., n)
# - plotting these statistics
# - doing t-tests between samples

import pdb

import numpy as np
import scipy.stats
import pandas

def event_stats(events, group_by):
    """Computes mean and s.e.m. of event statistics.

    events: pandas DataFrames of events from one experiment,
        as returned by get_events_at_region()
    group_by: the columns to group by
    Returns: a data frame including mean, s.d., and s.e.m. of
        current and dwell time.
    """
    events = events[ events.central_base != 'N' ]
    events_by_pos = events.groupby(group_by)
    event_summary = pandas.DataFrame({
        'current_mean': np.round(events_by_pos['event_level_mean'].mean(), 4),
        'current_std': np.round(events_by_pos['event_level_mean'].std(), 8),
        'current_sem': np.round(events_by_pos['event_level_mean'].sem(), 8),
        'dwell_time_mean': np.round(events_by_pos['dwell_time'].mean(), 5),
        'dwell_time_std': np.round(events_by_pos['dwell_time'].std(), 10),
        'dwell_time_sem': np.round(events_by_pos['dwell_time'].sem(), 10),
        'count': events_by_pos.size()
    })
    return event_summary

def event_stats_by_position(events):
    """Gets events by position."""
    return event_stats(events, ['contig', 'pos1', 'reference_kmer'])

def event_stats_by_kmer(events):
    """Gets events only by k-mer."""
    return event_stats(events, ['reference_kmer'])

def get_events(read_events, region):
    """Gets events from several read event data sets.

    read_events: hash with keys being sample names,
        and values being ReadEventData objects
    region: region at which to get events
    Returns: DataFrame of events, with an additional sample_name column
    """
    def get_events(sample_name):
        events = datasets[sample_name].get_events_at_region(region)
        events['sample_name'] = sample_name
        return events
    event_samples = [get_event_sample(sample_name)
        for sample_name in datasets.keys()]
    event_samples = pandas.concat(event_samples)

def t_test(mod_events, control_events):
    """Computes a t-test between two sets of events.

    mod_events: pandas DataFrames of events for modified sample,
        as returned by get_events_at_region()
    control_events: similarly, for control sample
    Returns: t-test, comparing current and dwell time between
        events from modified samples, versus unmodified.
    """
    # tag events by sample_name, and combine them
    mod_events['modified'] = True
    control_events['modified'] = False
    events = pandas.concat([mod_events, control_events])
    events = events[ events.central_base != 'N' ]
    t_test_results = []

    for stat_column in ['event_level_mean', 'dwell_time']:
        # convert events into a table, with one column for each position
        x = events.pivot(index=['modified', 'read_name'],
            columns=['contig', 'pos1'],
            values=stat_column)
        # do t-test
        r = scipy.stats.ttest_ind(
            x.loc[True].to_numpy(), x.loc[False].to_numpy(),
            equal_var=True, nan_policy='omit')
        positions = x.columns.to_frame()
        r = pandas.DataFrame({
            'stat': stat_column,
            'contig': positions.contig,
            'pos1': positions.pos1,
            't': r.statistic,
            'p': r.pvalue
        })
        r['stat'] = stat_column 
        t_test_results.append(r)
    return pandas.concat(t_test_results)

def event_stats_summary(events):
    """Gets statistics about some events.

    events: events, as returned by get_events_at_region(), including 'dwell_time'
    Returns: a data frame of stats about the event current and dwell time
    """
    events['current'] = events['event_level_mean']
    events['log10_dwell_time'] = np.log10(events['dwell_time'] * 1000.)
    # discard very-long events
    events = events[ events['dwell_time'] < 0.5]
    events = events[['current', 'log10_dwell_time']]
    # get statistics
    stats = events.describe(percentiles=[0.5])
    return stats

