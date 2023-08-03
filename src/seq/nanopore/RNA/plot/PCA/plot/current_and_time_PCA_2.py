#!/usr/bin/env python3
# Plots PCA of reads with various modifications.
# Also, should probably plot all the modifications.

import os
import pdb

import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
import pandas

from scipy.linalg import sqrtm

# read in PCA from Prism
# (this may change to using sklearn)
# pca = pandas.read_table('../PCA of Current and dwell time 1.txt')
# using sklearn
pca = pandas.read_csv('../PCA_columns_standardized.csv')
# add some annotations
pca['mod_name'] = [s.split(' ')[0] for s in pca['sample name']]
pca['concentration'] = [s.split(' ')[1] for s in pca['sample name']]

# get means for PCA from each method
PC_mean_prism = pandas.read_csv('../Prism/PC_scores_mean.csv')
PC_mean_prism.rename(columns={'Label': 'sample name'}, inplace=True)
PC_mean_scipy = pca.groupby('sample name')[['PC1','PC2','PC3']].mean()
PC_mean_scipy.reset_index(inplace=True)

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
        '100%': 5
        }
mod_names = list(set(pca['mod_name']))
mod_names.sort()
# set hues for the modifications
hue = dict(zip(mod_names, np.linspace(0, 1, len(mod_names), endpoint=False)))
# get color scheme for these
sample_color = dict()
for i in pca[['sample name', 'mod_name', 'concentration']].drop_duplicates().iterrows():
    sample_color[i[1]['sample name']] = matplotlib.colors.hsv_to_rgb((
            hue[i[1]['mod_name']],
            0.75,
            concentration_levels[i[1]['concentration']] / 5))

# get bounds (so that all the plots can be on the same scale)
pca_bounds = pca.loc[:,['PC1','PC2','PC3']].quantile([0,1])




def std_ellipse(sigma, num_points=1000):
    """Gets points corresponding to a 1-sd ellipse.

    Adapted from (Randolf Scholz' comment)[https://gist.github.com/CarstenSchelp/b992645537660bda692f218b562d0712?permalink_comment_id=3465086#gistcomment-3465086]
    sigma: the covariance matrix
    num_points: how many points in a circle to use
    Returns: a numpy array with shape (2, num_points), of
        the points on the ellipse.
    """
    T = np.linspace(0, 2*np.pi, num=num_points)
    circle = radius * np.vstack([np.cos(T), np.sin(T)])
    x, y = sqrtm(Sigma) @ circle
#    plt.plot(x, y, '-r', linewidth=5)


def plot_center_and_conf_ellipse(mu, sigma, scale=2.44, color='black', lw=2):
    """

    mu, sigma: mean and covariance
    scale: the scale of the ellipse (relative to s.d. along each major axis)
        According to the link below, setting the half-axis scale to
        2.44 = sqrt(5.991) (approximately) will give an ellipse
        covering 95% of the points.
        https://www1.udel.edu/biology/rosewc/kaap686/reserve/cop/center%20of%20position%20conf95.pdf
    color: color to use
    r: mean of radius
    lw: width of line to use
    Side effects: plots the mean (as a point) and the (scaled) s.d. ellipse
    """
    plt.scatter(mu, c=color, r=lw)
    plt.plot(mu + scale * std_ellipse(sigma), c=color, lw=lw)


def plot_PCA_for_mod(mod_name):
    """Plots PCA of all the concentrations of one modification.

    mod_name: the name of the modification to plot
    Side effects: plots PCA for that modification
    """
    hue = 2/3   # FIXME set this in advance
    # using the first two components seems to work reasonably well
    components = ['PC1', 'PC2']
    os.makedirs(output_dir, exist_ok=True)
    # set up plotting
    plt.figure(figsize=(9,6))
    for concentration in concentration_levels.keys():
        x1 = pca[(pca.mod_name==mod_name) & (pca.concentration==concentration)]
        if x1.empty:
            continue
        # pdb.set_trace()
        sample_name = f'{mod_name} {concentration}'
        plt.scatter(x1[components[0]], x1[components[1]],
                label=sample_name,
                color=sample_color[sample_name],
                alpha=0.8)
        # FIXME plot just the centers of these
    # label axes
    plt.xlabel(components[0])
    plt.ylabel(components[1])
    # set bounds (to be consistent across plots)
    xlim=pca_bounds[components[0]]
    plt.xlim(xlim[0], xlim[1])
    ylim=pca_bounds[components[1]]
    plt.ylim(ylim[0], ylim[1])
    # FIXME add legend?
    plt.legend()
    plt.savefig(f'{output_dir}/{mod_name}.png')

def plot_all_sample_mean(PC_means, output_prefix):
    """Plots mean for all of the samples."""
    for x_label, y_label in [('PC1','PC2'), ('PC2','PC3')]:
        plt.figure(figsize=(7,7))
        plt.scatter(PC_means[x_label], PC_means[y_label],
                color='black', s=5, alpha=0.7)
        for i in range(len(PC_means)):
            plt.text(PC_means.loc[i,x_label],
                   PC_means.loc[i,y_label],
                   ' ' + PC_means.loc[i,'sample name'],
                   color='black', alpha=0.7)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.savefig(f'{output_prefix}_{x_label}_{y_label}.png')

if False:
    for mod in mod_names:
        print(mod)
        plot_PCA_for_mod(mod)

plot_all_sample_mean(PC_mean_prism, 'PCA_prism')
plot_all_sample_mean(PC_mean_scipy, 'PCA_scipy')

