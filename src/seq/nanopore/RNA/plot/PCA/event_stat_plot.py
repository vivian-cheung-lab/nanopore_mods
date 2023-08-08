#!/usr/bin/env python3
# Plots distribution of f-statistics,
# and stats about columns with high f-statistics.

import pdb
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas
import seaborn as sns

# import pysam

stats = pandas.read_csv('read_event_stats_ANOVA_bases_5000reads.csv.gz')
stats.set_index(stats.columns[0], inplace=True)
# stats = stats.iloc[:,1:]

stats['measurement'] = stats.index.str.replace(' \d+', '', regex=True)

# extract stats (in this case, f-statistic from one-way ANOVA)
f_stats = pandas.DataFrame({
    'Current f-stat': stats[stats.measurement=='current'].statistic.tolist(),
    'Dwell time f-stat': stats[stats.measurement=='time'].statistic.tolist()
    })

quantile_cutoff = 0.05
current_cutoff = np.quantile(f_stats.iloc[:,0], 1-quantile_cutoff)
time_cutoff = np.quantile(f_stats.iloc[:,1], 1-quantile_cutoff)

# plot distribution of f-statistics
plt.figure(figsize=(10,4))
fig, axes = plt.subplots(nrows=1, ncols=2)

# pdb.set_trace()

z = pandas.DataFrame({'x': [1,2,3,4,5]})
# N.B. this seems to require seaborn >= 0.12.2
sns.histplot(f_stats.iloc[:,0], ax=axes[0], bins=20)
axes[0].axvline(x=current_cutoff, color='black')
sns.histplot(f_stats.iloc[:,1], ax=axes[1], bins=20)
axes[1].axvline(x=time_cutoff, color='black')
plt.tight_layout()
plt.savefig('event_stat_plot.pdf')

##### plot info about selected columns

print('[reading in events]')
events = pandas.read_csv('read_event_table_bases_5000reads.csv.gz', nrows=20000)  # for practice
events_by_sample = events.set_index('sample name').iloc[:,1:]
sample_names = events_by_sample.index.unique()

def plot_distributions(output_pdf_name, column_stats):
    """Plots distribution of a column.

    Possibly not used.
    output_pdf_filename: where to write figures
    column_stats: info about the columns to include
    Returns: a figure with the distribution of that
    """
    pdf = matplotlib.backends.backend_pdf.PdfPages(output_pdf_name)
    for i in column_stats.itertuples():
        fig = plt.figure()
        pdb.set_trace()
        x = events_by_sample[i]
        # sns.histplot(data=x, hue=x.index,
        # element='step', fill=False)

        pdf.savefig(fig)
    pdf.close()



events1 = events_by_sample      # .iloc[:20000,:5]

stats = events1.groupby('sample name').describe(percentiles=[0.5])
# only include first copy of "counts" column
counts = stats.iloc[:,0]
counts.columns = ['n']
stats = pandas.concat([counts, stats.loc[:, [i[1]!='count' for i in stats.columns]] ], axis=1) 
stats = stats.T

stats.to_excel('read_event_stats.xlsx')
# plot_distributions('read_event_column_distributions.pdf', stats[:5])


# get clone sequence
# ref_seq = pysam.FastaFile(ref_fasta_file)




