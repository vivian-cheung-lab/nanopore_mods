#!/usr/bin/env python3

import pdb

import matplotlib.pyplot as plt
import numpy as np
import pandas
import seaborn as sns

stats = pandas.read_csv('read_event_stats_ANOVA_bases_5000reads.csv.gz')
stats.set_index(stats.columns[0], inplace=True)
# stats = stats.iloc[:,1:]

stats['measurement'] = stats.index.str.replace(' \d+', '')

# extract stats (in this case, f-statistic from one-way ANOVA)
f_stats = pandas.DataFrame({
    'Current f-stat': stats[stats.measurement=='current'].statistic.tolist(),
    'Dwell time f-stat': stats[stats.measurement=='time'].statistic.tolist()
    })

quantile_cutoff = 0.05
current_cutoff = np.quantile(f_stats.iloc[:,0], 1-quantile_cutoff)
time_cutoff = np.quantile(f_stats.iloc[:,1], 1-quantile_cutoff)


# plot distribution of standard deviations
plt.figure(figsize=(10,4))
fig, axes = plt.subplots(nrows=1, ncols=2)
sns.histplot(f_stats.iloc[:,0], ax=axes[0], bins=20)
axes[0].axvline(x=current_cutoff, color='black')
sns.histplot(f_stats.iloc[:,1], ax=axes[1], bins=20)
axes[1].axvline(x=time_cutoff, color='black')
plt.tight_layout()
plt.savefig('event_stat_plot.pdf')

