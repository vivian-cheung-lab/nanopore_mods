#!/usr/bin/env python
# coding: utf-8

# # Plotting error rates
# These are attempts at different ways of plotting error rates.
# (Adapted from a Jupyter notebook.)

# In[1]:


import re

import pandas

import pdb
import matplotlib.pyplot as plt
import seaborn as sns


# In[2]:


# location of data
data_base = '../'

# sample info
sample_table = pandas.read_csv(data_base + '/IVT_samples.csv')
sample_table.rename(columns={'FASTQ_file_base_name': 'Sample name'},
                    inplace=True)
# slight renaming
sample_table.short_name[sample_table.short_name=='biotin-C'] = 'Bio-C'
sample_table.short_name[sample_table.short_name=='Y'] = '𝚿'
# also tack on "label" (used for graphing)
sample_table['label'] = sample_table.mod_concentration + ' ' + sample_table.short_name

sample_table.set_index('Sample name', inplace=True)
sample_table


# ## Mismatch stats
# First, we look at mismatch stats.

# In[3]:


def get_mismatch_stats_summary(sample_name):
    """Gets summary of statistics about mismatches."""
    # open a dataset
    RDD = pandas.read_csv(data_base + '/RDD_3/' + sample_name + '.csv.gz')
    RDD.rename(columns = {'Ref': 'Seq'}, inplace=True)
    # count bases of each type
    counts = RDD.groupby('Seq')[['A','C','G','U','n']].sum()
    counts['Bases'] = counts[['A','C','G','U']].sum(axis=1)
    # count mismatches
    counts['Match'] = [counts.loc[b,b] for b in 'ACGU']
    counts['Mismatch'] = counts['Bases'] - counts['Match']
    # mismatches (relative to total number of bases; it's not
    # entirely clear what the denominator should be here, 'n' or 'Bases')
    counts['Mismatch_percent'] = 100 * counts['Mismatch'] / counts['n']
    
    stats = pandas.DataFrame(counts['Mismatch_percent']).transpose()
    stats['Sample name'] = sample_name
    # stats.index = [sample_name]
    # stats.columns.name = 'Sample name'
#    stats.reset_index(drop=True, inplace=True)
    return stats

mismatch_stats = pandas.concat([
    get_mismatch_stats_summary(s)
    for s in sample_table.index.to_list()])
mismatch_stats.columns.name = None
mismatch_stats.set_index('Sample name', inplace=True)


# In[4]:


mismatch_stats.columns.name = None
mismatch_stats.columns
mismatch_stats = sample_table.join(mismatch_stats)
mismatch_stats.head()


# Here are the mismatch stats for just the canonical samples.

# In[5]:


mismatch_stats_canonical = mismatch_stats[mismatch_stats.mod_concentration=='none']
mismatch_stats_canonical


# ## Mismatch stats for all N>X types

# We also break down the mismatch stats for each sort of substitution.

# In[6]:


def get_mismatch_stats_by_type(sample_name):
    """Gets summary of statistics for all X>Y types of mismatches."""
    # open a dataset
    RDD = pandas.read_csv(data_base + '/RDD_3/' + sample_name + '.csv.gz')
    RDD.rename(columns = {'Ref': 'Seq'}, inplace=True)
    # count bases of each type
    counts = RDD.groupby('Seq')[['A','C','G','U','n']].sum()
    # counts['Bases'] = counts[['A','C','G','U']].sum(axis=1)
    
    mismatch_percentage = {b: (100 * counts[b] / counts.n)
                           for b in 'ACGU'}
    mismatch_percentage = pandas.DataFrame(mismatch_percentage)
    mismatch_percentage = mismatch_percentage.unstack().transpose()
    mismatch_percentage.index = [f'{b[1]}>{b[0]}' for b in mismatch_percentage.index]
    mismatch_percentage = pandas.DataFrame(mismatch_percentage,
                                          columns=[sample_name])
    mismatch_percentage = mismatch_percentage.transpose()
    return mismatch_percentage

mismatch_stats_by_type = pandas.concat([
    get_mismatch_stats_by_type(s)
    for s in sample_table.index.to_list()])
mismatch_stats_by_type.index.name = 'Sample name'
mismatch_stats_by_type = sample_table.join(mismatch_stats_by_type)

mismatch_stats_by_type.head()


# Next, we get just the U>X stats, for the canonical samples.
# 

# In[7]:


canonical_U_to_X = mismatch_stats_by_type.loc[
    mismatch_stats_by_type.mod_concentration=='none',
    ['short_name', 'mod_concentration', 'label', 'U>A', 'U>C', 'U>G']]
canonical_U_to_X


# ### Mismatch stats for canonical samples
# 

# We'll need averages for the canonical bases.
# 

# In[8]:


mismatch_stats_canonical_mean = mismatch_stats_canonical[['A','C','G','U']].mean().to_frame().transpose()
# XXX this is a bit hokey
mismatch_stats_canonical_mean = pandas.concat(
    [pandas.DataFrame([{'short_name':'canonical', 'mod_concentration':'none', 'label':'canonical'}]),
    mismatch_stats_canonical_mean], axis=1)
mismatch_stats_canonical_mean.index = ['canonical']
mismatch_stats_canonical_mean


# In[9]:


def sorted_order(x):
    """Utility for sorting by various columns.
    
    Note that this includes various sorts of data.
    """
    ordering = ['none', '1 to 100', '1 to 10', '1 to 5', '1 to 2', '100%',
               'm1A', 'm6A', 'Bio-C', 'hm5C', 'm5C', '𝚿', 'Am', 'Cm', 'Gm', 'Um']
    d = dict(zip(ordering, range(len(ordering))))
    return [d[i] for i in x]


# ## Mismatch stats tables
# We format these to include the canonical stats.
# 
# First, the base stats.

# In[10]:


mismatch_stats_base = mismatch_stats[
    (~mismatch_stats.short_name.str.match('[ACGU]m'))
    & (~mismatch_stats.mod_concentration.str.match('(none|100%)'))].copy()
# pdb.set_trace()
mismatch_stats_base.sort_values(['short_name', 'mod_concentration'], key=sorted_order, inplace=True)
# the above, with canonicals prepended
mismatch_stats_base_1 = pandas.concat([
    mismatch_stats_canonical_mean, mismatch_stats_base
])
mismatch_stats_base_1


# In[11]:


mismatch_stats_sugar = mismatch_stats[
    (mismatch_stats.short_name.str.match('[ACGU]m'))
    & (~mismatch_stats.mod_concentration.str.match('(none)'))].copy()
mismatch_stats_sugar.sort_values(['short_name', 'mod_concentration'], key=sorted_order, inplace=True)
# the above, with canonicals prepended
mismatch_stats_sugar_1 = pandas.concat([
    mismatch_stats_canonical_mean, mismatch_stats_sugar
])
mismatch_stats_sugar_1


# In[12]:


mismatch_stats_sugar_100pct = mismatch_stats[
    (mismatch_stats.short_name.str.match('[ACGU]m'))
    & (mismatch_stats.mod_concentration=='100%')].copy()
mismatch_stats_sugar_100pct.sort_values(['short_name', 'mod_concentration'], key=sorted_order, inplace=True)
# the above, with canonicals prepended
mismatch_stats_sugar_100pct = pandas.concat([
    mismatch_stats_canonical_mean, mismatch_stats_sugar_100pct
])
mismatch_stats_sugar_100pct


# ## Indel stats
# Next, we look at indel stats.

# In[13]:


# get indel stats
def get_indel_stats_summary(sample_name):
    """Gets summary of statistics about indels."""
    # open a dataset
    RDD = pandas.read_csv(data_base + '/RDD_3/' + sample_name + '.csv.gz')
    # RDD['Sample'] = RDD['Sample'].str.replace('\.bam$', '', regex=True)
    # return RDD
    # RDD.set_index('Sample', inplace=True)
    n = RDD.n.sum()
    return pandas.DataFrame([{
        'Sample name': sample_name,
        'Insertion sites': 100 * RDD.I.sum() / n,
        'Deletion sites': 100 * RDD.X.sum() / n,
        'Inserted bases': 100 * RDD.i.sum() / n,
        'Deleted bases': 100 * RDD.D.sum() / n
    }])

indel_stats = pandas.concat([
    get_indel_stats_summary(s)
    for s in sample_table.index.to_list()])
indel_stats.set_index('Sample name', inplace=True)
indel_stats = indel_stats[['Inserted bases', 'Deleted bases']]
indel_stats = sample_table.join(indel_stats)
indel_stats.head()


# ### Indel stats for canonical samples
# As for the mismatches, we'll need the average indels for the canonical samples.

# In[14]:


indel_stats_canonical = indel_stats[indel_stats.mod_concentration=='none']
indel_stats_canonical


# In[15]:


indel_stats_canonical_mean = indel_stats_canonical[['Inserted bases', 'Deleted bases']].mean().to_frame().transpose()
# XXX this is a bit hokey
indel_stats_canonical_mean = pandas.concat(
    [pandas.DataFrame([{'short_name':'canonical', 'mod_concentration':'none', 'label':'canonical'}]),
    indel_stats_canonical_mean], axis=1)
indel_stats_canonical_mean.index = ['canonical']
indel_stats_canonical_mean


# ### Indel stats tables
# 
# Here are the indel error rates for the base modifications (not including 100%).

# In[16]:


indel_stats_base = indel_stats[
    (~indel_stats.short_name.str.match('[ACGU]m'))
    & (~indel_stats.mod_concentration.str.match('(none|100%)'))].copy()
indel_stats_base.sort_values(['short_name', 'mod_concentration'], key=sorted_order, inplace=True)
# the above, with canonicals prepended
indel_stats_base_1 = pandas.concat([
    indel_stats_canonical_mean, indel_stats_base
])
indel_stats_base_1


# Here are the indel error rates for the sugar modifications (not including 100%).

# In[17]:


indel_stats_sugar = indel_stats[
    indel_stats.short_name.str.match('[ACGU]m')
    & ~indel_stats.mod_concentration.str.match('(none)')].copy()
indel_stats_sugar.sort_values(['short_name', 'mod_concentration'], key=sorted_order, inplace=True)
# the above, with canonicals prepended
indel_stats_sugar_1 = pandas.concat([
    indel_stats_canonical_mean, indel_stats_sugar
])
indel_stats_sugar_1.head()


# ... and the 100% modified.

# In[18]:


indel_stats_sugar_100pct = indel_stats[
    indel_stats.short_name.str.match('[ACGU]m')
    & (indel_stats.mod_concentration=='100%')].copy()
indel_stats_sugar_100pct.sort_values(['short_name', 'mod_concentration'], key=sorted_order, inplace=True)
# the above, with canonicals prepended
indel_stats_sugar_100pct = pandas.concat([
    indel_stats_canonical_mean, indel_stats_sugar_100pct
])
indel_stats_sugar_100pct.head()


# # Sheets in Excel format

# In[19]:

# for now, omitting the "100% sugar" tables
with pandas.ExcelWriter('error_rates_2.xlsx') as xls:
    mismatch_stats_canonical.to_excel(xls, sheet_name='Mismatch, canonical')
    canonical_U_to_X.to_excel(xls, sheet_name='U>X, canonical')
    indel_stats_canonical.to_excel(xls, sheet_name='Indels, canonical')
    
    mismatch_stats_base_1.to_excel(xls, sheet_name='Mismatch, bases (<100%)')
    mismatch_stats_sugar_1.to_excel(xls, sheet_name='Mismatch, sugars')
    # mismatch_stats_sugar_100pct.to_excel(xls, sheet_name='Mismatch, sugars (100%)')
    
    indel_stats_base_1.to_excel(xls, sheet_name='Indels, bases (<100%)')
    indel_stats_sugar_1.to_excel(xls, sheet_name='Indels, sugars')
    # indel_stats_sugar_100pct.to_excel(xls, sheet_name='Indels, sugars (100%)')


# In[ ]:




