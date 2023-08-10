#!/usr/bin/env python3
# Reformats the ANOVA results slightly, and tacks on the central base.

import pdb

import numpy as np
import pandas

import pysam

# get clone sequence
clone_fasta = pysam.FastaFile('../../../../../../chr19_AANCR_clone_only_230411.fa')
# note that this is AANCR's sequence, with 0-based coordinates
clone_seq = clone_fasta.fetch('chr19', 45406984, 45408892).upper()

ANOVA_stats = pandas.read_csv('read_event_stats_ANOVA_bases_5000reads.csv.gz')
ANOVA_stats['measurement'] = [s[0]
        for s in ANOVA_stats.iloc[:,0].str.split()]
ANOVA_stats['base_number'] = [int(s[1])
        for s in ANOVA_stats.iloc[:,0].str.split()]

ANOVA_table = ANOVA_stats.pivot(index='base_number', columns='measurement')
ANOVA_table = ANOVA_table.iloc[:,2:]
# note that the base index in the table is 1-based
ANOVA_table['base'] = [clone_seq[(i-3):(i+2)].replace('T', 'U')
        for i in ANOVA_table.index]

def base_composition_at_cutoff(proportion_to_include):
    """Gets base composition of the sites with the top % of F-statistics."""
    def above_cutoff(x):
        cutoff = np.quantile(x, 1-proportion_to_include)
        return x >= cutoff
    bases = ANOVA_table[
            above_cutoff(ANOVA_table.statistic.current) |
            above_cutoff(ANOVA_table.statistic.time) ].base.str[2]
    base_composition = bases.to_frame().groupby('base').size().to_frame()
    base_composition['cutoff'] = proportion_to_include
    return base_composition

base_composition = pandas.concat([base_composition_at_cutoff(cutoff)
    for cutoff in [0.01, 0.05, 0.1, 0.25, 0.5]])
# pdb.set_trace()

base_composition.reset_index(inplace=True)
base_composition = base_composition.pivot(index=['cutoff'], columns='base')
base_composition.columns = base_composition.columns.droplevel(0)
# pdb.set_trace()

with pandas.ExcelWriter('event ANOVA and base.xlsx', engine='xlsxwriter') as writer:
    ANOVA_table.to_excel(writer, sheet_name='ANOVA and central base')
    base_composition.to_excel(writer, sheet_name='Base composition')

