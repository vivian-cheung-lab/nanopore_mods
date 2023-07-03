#!/usr/bin/env python3
# Plots event statistics at particular bases.
# FIXME
# - this is not that fast
# - change font of '𝚿' to something besides Arial?
#   (possibly not worth the effort)

import itertools
import multiprocessing
import os
import pdb
import subprocess
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas
import pysam
import seaborn as sns

# add local modules to this
sys.path.append('py/')

import plot.seaborn.utils
from statannotations.Annotator import Annotator

import seq.nanopore.squiggle.read_event

# set default font to Arial
plt.rcParams['font.family'] = 'Arial'

# where MOP2 output is
MOP2_dir = '../data/AANCR_IVT_1/mop_preprocess/output/'

# root of the other files
data_base = './'

# where to write plots
output_dir = 'plots/'

os.makedirs(output_dir, exist_ok=True)

# where to write out events to (omitting for now)
event_tables_dir = None     # 'event_tables/'

# where to write event stats to
event_stats_dir = 'output/'

os.makedirs(output_dir, exist_ok=True)
os.makedirs(event_stats_dir, exist_ok=True)

def get_dataset(sample_name):
    """Gets a dataset."""
    # open a dataset 
    read_event = seq.nanopore.squiggle.read_event.ReadEventData(
        MOP2_dir + '/fastq_files/' + sample_name + '.fq.gz',
        data_base + '/events/' + sample_name + '.mergedEvent.tsv.gz')
    return read_event

# region = 'chr19:45406985-45408892'
# region = 'chr19:45407000-45407300'   # for practice
# practice_region = 'chr19:45407000-45407050'   # for practice
# for 100% sugar modifications, only including this region
# (as that's around where their read depth drops off)
# region = 'chr19:45406985-45407485'

sample_table = pandas.read_csv(data_base + '/IVT_samples.csv')
# sample_table.rename(columns={'FASTQ_file_base_name': 'Sample name'},
#                     inplace=True)
# sample_table.set_index('Sample name', inplace=True)

def get_sample_name(short_name, mod_concentration):
    """Gets sample name by short name and concentration.

    short_name, mod_concentration: these specify which sample to get
    Returns: the relevant sample name, or None on failure
    """
    print(f'{short_name} {mod_concentration}')
    sample_name = sample_table[(sample_table.short_name==short_name)
        & (sample_table.mod_concentration==mod_concentration)].FASTQ_file_base_name
    if sample_name.shape[0] == 1:
        return sample_name.iloc[0]
    return None

def get_events_by_sample_name(sample_name, region):
    """Gets events for one sample."""
    fastq = f'{data_base}/fastq/{sample_name}.fastq.gz'
    events_file = f'{data_base}/events/{sample_name}.mergedEvent.tsv.gz'
    read_data = seq.nanopore.squiggle.read_event.ReadEventData(
        fastq, events_file)
    events = read_data.get_events_at_region(region,
        include_dwell_time=True)
    return events

def get_events(short_name, mod_concentration, region):
    """Gets events for one pair of samples."""
    mod_events = get_events_by_sample_name(
        get_sample_name(short_name, mod_concentration), region)
    mod_events['Sample'] = 'Modified'
    mod_events['Modification'] = short_name
    mod_events['Concentration'] = mod_concentration
    unmod_events = get_events_by_sample_name(
        get_sample_name(short_name, 'none'), region)
    unmod_events['Sample'] = 'Unmodified'
    unmod_events['Modification'] = short_name
    unmod_events['Concentration'] = 'none'
    events = pandas.concat([mod_events, unmod_events])
    # XXX hack to set name for pseudouridine
    if short_name == 'Y':
        events['Modification'] = '𝚿'
    # tack on dwell time in ms
    events['dwell_time_ms'] = 1000. * events['dwell_time']
    # note that this is log10( dwell time in ms )
    events['log10_dwell_time'] = np.log10(events['dwell_time'] * 1000.)
    return events

# XXX custom p-value thresholds (for use with "star" style of p-value formatting)
custom_p_value_thresholds = (
    [[10 ** (-i), f'p<1e-{i}'] for i in range(300, 2, -1)]
    + [
    [1e-2, 'p<0.01'],
    [0.05, 'p<0.05'],
    [1., '']])

def get_data(short_name, mod_concentration, region, stat_name, center_of_kmer):
    """Gets event data.

    short_name: short name of the modification
    mod_concentration: which concentration to use
    region: the region at which to get events
    stat_name: either 'Current', 'Dwell time', or 'Log10(Dwell time)'
    center_of_kmer: if True, restrict to events with the base in
        the center of the kmer. If False, restrict to events with
        the given base _anywhere_ in the kmer.
    Returns: data about that stat
    """
    events = get_events(short_name, mod_concentration, region)
    # discard outliers
    events = events[ events['dwell_time'] < 0.5]
    # for now, getting stats for all central bases (so as to compare them),
    # even though we probably won't include all of them
    # pdb.set_trace()
    def stats_at_base(b):
        if center_of_kmer:
            events1 = pandas.DataFrame(events[events.central_base==b])
        else:
            events1 = pandas.DataFrame(events[events.reference_kmer.str.contains(b)])
        events1['Central base'] = b
        print(events1.shape)
        return events1
    events_by_base = pandas.concat([stats_at_base(b) for b in 'ACGT'])
    events_by_base['Central base'] = events_by_base['Central base'].str.replace('T', 'U')
    return events_by_base

def plot_dist_comparison(events_by_base, stat_name, plot_width, output_name, mod_color='red'):
    """Plots a comparison, for each different base.
   
    """
    # XXX hack to set name for pseudouridine
    # note that this is in Arial font; using a different font seems difficult
    events_by_base['Modification'].replace('𝚿', '$\Psi$', inplace=True)
    center_of_kmer = True
    column_name = {'Current': 'event_level_mean',
        'Dwell time': 'dwell_time',
        'Dwell time (ms)': 'dwell_time_ms',
        'Log10(Dwell time)': 'log10_dwell_time'}[stat_name]
    modifications = list(set(events_by_base['Modification']))
    modifications.sort()
    fig = plt.figure(figsize=(plot_width, 4.8))
    # style, for setting lines to black
    s = {'color': 'black'}
    # plot parameters (split out for using statsannotation)
    plot_params = {
        'data': events_by_base,
        'x': 'Modification',
        'y': column_name,
        'hue': 'Sample',
        'sym': '',            # omit plotting outliers
        'hue_order': ['Unmodified', 'Modified'],

        # workaround for issue with seaborn using np.float
        # (which numpy has deprecated)
        'orient': 'v',
        'order': modifications,

        # color bars
        'palette': {'Unmodified': 'grey', 'Modified': mod_color},
        # set lines in these to black
        'capprops': s, 'whiskerprops': s, 'medianprops': s,
        'boxprops': {'edgecolor': 'black'}
    }
# FIXME
# - tweak y-axis?
    bp = sns.boxplot(**plot_params)
    # plt.title(output_name.replace('_', ' '))
    # plt.xlabel('Base at center of 5-mer' if center_of_kmer else 'Base anywhere in 5-mer')
    plt.xlabel('Modification')
    # plt.ylabel('Dwell time (s)')
    # omitting title
    # plt.title(output_name)
    label = 'base_at_center_of_kmer' if center_of_kmer else 'base_in_kmer'
    # add space between boxes in the same column
    plot.seaborn.utils.adjust_box_widths(fig, 0.75)
    # this is based on example code at
    # https://levelup.gitconnected.com/statistics-on-seaborn-plots-with-statannotations-2bfce0394c00
    pairs = tuple([[(m, 'Unmodified'), (m, 'Modified')] for m in modifications])
    # by default, using human-readable name for this
    plt.ylabel(stat_name)
    if stat_name == 'Dwell time (ms)':
        plt.ylabel('Dwell time (milliseconds)')
    # possibly set y-axis limits
    if stat_name == 'Log10(Dwell time)':
        # XXX hack to label with log10-scale
        plt.ylabel('Dwell time (ms, log 10)')
        bp.set_ylim(0, 2.5)
        bp.yaxis.set_ticks([0, 1, 2], labels=['1', '10', '100'])
    # add annotations
    annotator = Annotator(bp, pairs, **plot_params)
    annotator.configure(test='t-test_ind', comparisons_correction='bonferroni',
        text_format='star', show_test_name=False)
#        pvalue_thresholds=custom_p_value_thresholds)
        # avoid printing *'s for non-significant cases?
        # pvalue_thresholds=[[1e-4, "****"], [1e-3, "***"], [1e-2, "**"], [0.05, "*"], [1, ""]])
        # [[1e-4, "****"], [1e-3, "***"], [1e-2, "**"], [0.05, "*"], [1, ""]])
    # for now, just printing stars
    # ??? show p-values?
    _, corrected_results = annotator.apply_and_annotate()
    print(corrected_results)

    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{output_name}.pdf')
    return corrected_results

def plot_comparison(output_name, concentration, region, plot_width, mods_and_bases, mod_color='red'):
    """Plots a comparison of some modifications.

    """
    def get_mod_events(mod_and_base):
        modification, base = mod_and_base
        events = get_data(modification, concentration, region, 'Dwell time', True)
        # restrict to only the events with that modification (possibly) at the center
        events = events[events['Central base']==base]
        return events
    events = pandas.concat([get_mod_events(mb) for mb in mods_and_bases])
    # possibly write out tables of events
    if event_tables_dir:
        os.makedirs(event_tables_dir, exist_ok=True)
        events.to_csv(f'{event_tables_dir}/Events {output_name}.csv.gz')
    # possibly write out event stats
    if event_stats_dir:
        os.makedirs(event_stats_dir, exist_ok=True)
        # events['Dwell time'] = events.dwell_time
        events['Log10(Dwell time)'] = events['log10_dwell_time']
        # omitting log10-transformed stats
        if False:
            event_stats = events.groupby(
                ['Modification', 'Concentration', 'Sample'])['Log10(Dwell time)'].describe(percentiles=[0.5])
            event_stats = event_stats.round(5)
            event_stats.rename(columns={'50%': 'median'}, inplace=True)
            event_stats['count'] = event_stats['count'].astype(int)
        # using non-log-transformed stats
        event_stats = events.groupby(
            ['Modification', 'Concentration', 'Sample'])['dwell_time_ms'].describe(percentiles=[0.5])
        event_stats = event_stats.round(1)
        event_stats.rename(columns={'50%': 'median'}, inplace=True)
        event_stats['count'] = event_stats['count'].astype(int)
        event_stats.to_csv(f'{event_stats_dir}/Dwell time for {concentration} {output_name}.csv.gz')
    plot_dist_comparison(events, 'Dwell time (ms)', plot_width,
        f'Dwell time for {output_name}, {concentration}', mod_color=mod_color)

if __name__ == '__main__':
    # plot_sugar_comparison('1:2')
    # pool = multiprocessing.Pool(30)
    # stats = pool.map(plot_sugar_comparison, ['100%', '1:2', '1:5', '1:10', '1:100'])
    # ??? this avoids an error message about NoneType
    # pool.close()

    # possibly plot only a few, for testing
    if False:
        plot_comparison('base modifications', '1 to 2', 'AANCR_IVT:1-1908', 7.9,
            [
            ('m1A', 'A')
            # ('biotin-C', 'C'),
            # ('Y', 'U'),
            ],
            mod_color='green')
        sys.exit(0)
    if True:
        plot_comparison('base modifications', '1 to 2', 'AANCR_IVT:1-1908', 7.9,
            [
            ('m1A', 'A'),
            ('m6A', 'A'),
            ('biotin-C', 'C'),
            ('hm5C', 'C'),
            ('m5C', 'C'),
            ('Y', 'U'),
            ],
            mod_color='green')
        plot_comparison('sugar modifications', '100%', 'AANCR_IVT:1-500', 6.4,
            [('Am', 'A'), ('Cm', 'C'), ('Gm', 'G'), ('Um', 'U')],
            mod_color='red')
    if False:
      plot_comparison('sugar modifications', '1:2', 'chr19:45406985-45407485', 6.4,
        [('Am', 'A'), ('Cm', 'C'), ('Gm', 'G'), ('Um', 'U')],
        mod_color='orange')

