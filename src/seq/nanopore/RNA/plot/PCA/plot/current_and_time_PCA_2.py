#!/usr/bin/env python3
# Plots PCA of reads with various modifications.

import os
import pdb

import matplotlib.pyplot as plt
import pandas

# read in PCA from Prism
# (this may change to using sklearn)
# pca = pandas.read_table('../PCA of Current and dwell time 1.txt')
# using sklearn
pca = pandas.read_csv('../PCA_columns_standardized.csv')

# where to write output
output_dir = 'current_and_time_PCA_2'

# mapping from concentration to "level of modification",
# for color-coding
concentration_levels = {
        'none': 0,
        '1:100': 1,
        '1:10': 2,
        '1:5': 3,
        '1:2': 4,
#        '100%': 5
        }



# get bounds (so that all the plots can be on the same scale)
pca_bounds = pca.loc[:,['PC1','PC2','PC3']].quantile([0,1])

mod_names = list(set([s.split(' ')[0] for s in set(pca['sample name'])]))
mod_names.sort()

def plot_PCA_for_mod(mod_name, hue):
    """Plots PCA of all the concentrations of one modification.

    mod_name: the name of the modification
    hue: the hue to use for plotting this modification
    Side effects: plots PCA for that modification
    """
    # using the first two components seems to work reasonably well
    components = ['PC1', 'PC2']
    os.makedirs(output_dir, exist_ok=True)
    x = pca.loc[:,['sample name',components[0],components[1]]]
    x = x[ x['sample name'].str.startswith(mod_name) ]
    x['concentration'] = x['sample name'].str.replace(mod_name + ' ', '')
    # set up plotting
    plt.figure(figsize=(9,6))
    for concentration in concentration_levels.keys():
        x1 = x[x.concentration==concentration]
        plt.scatter(x1[components[0]], x1[components[1]],
                label=concentration, alpha=0.4)
    plt.xlabel(components[0])
    plt.ylabel(components[1])
    plt.savefig(f'{output_dir}/{mod_name}.png')

for mod in mod_names:
    print(mod)
    plot_PCA_for_mod(mod, 2/3)

