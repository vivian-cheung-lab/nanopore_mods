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

# fraction of most-significant columns to include
quantile_cutoff = 0.05

def above_cutoff(x):
    """Gets only the rows with 'statistic' column above a cutoff."""
    cutoff = np.quantile(x.statistic, 1-quantile_cutoff)
    return x[ x.statistic >= cutoff ]

ANOVA_table = pandas.concat([
    above_cutoff(ANOVA_stats[ ANOVA_stats.measurement=='current' ]),
    above_cutoff(ANOVA_stats[ ANOVA_stats.measurement=='time' ]) ])
# note that the base index in the table is 1-based
ANOVA_table['bases'] = [clone_seq[(i-3):(i+2)].replace('T', 'U')
        for i in ANOVA_table.base_number]

print('base composition:')
print(ANOVA_table.bases.str[2].to_frame().groupby('bases').size())

with pandas.ExcelWriter('event ANOVA.xlsx', engine='xlsxwriter') as writer:
    ANOVA_table.to_excel(writer, sheet_name='Column f-statistics')

