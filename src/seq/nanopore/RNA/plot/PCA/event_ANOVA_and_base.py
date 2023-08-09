#!/usr/bin/env python3
# Reformats the ANOVA results slightly, and tacks on the central base.

import pdb

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
ANOVA_table['base'] = [clone_seq[i-1] for i in ANOVA_table.index]

ANOVA_table.to_xlsx('event_ANOVA_and_base.xlsx')

