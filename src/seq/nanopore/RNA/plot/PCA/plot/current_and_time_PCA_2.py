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
# omit sugars
# FIXME this is essentially a roundabout way of getting the sample names
pca = pca[ ~ pca.mod_name.str.match('^[A-Z]m$') ]

# where to write output
output_dir = 'current_and_time_PCA_2_bases_5000reads'

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


def plot_PCA_for_mod(pca, output_prefix, mod_name):
    """Plots PCA of all the concentrations of one modification.

    FIXME add all the points?
    pca: DataFrame of PCA scores
    output_prefix: output prefix for this
    mod_name: the name of the modification to plot
        (or None to include all of the modifications)
    Side effects: plots PCA for that modification
    """
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    # get bounds (so that all the plots can be on the same scale)
    pca_bounds = pca.loc[:,['PC1','PC2','PC3']].quantile([0,1])
    # using the first two components seems to work reasonably well
    components = ['PC1', 'PC2']
    # set up plotting
    plt.figure(figsize=(9,6))
    # include either this modification (or all of them)
    for modification in [mod_name] if mod_name else mod_names:
        for concentration in concentration_levels.keys():
            print(f'modification = {modification}')
            x1 = pca[(pca.mod_name==modification) & (pca.concentration==concentration)]
            if x1.empty:
                continue
            sample_name = f'{modification} {concentration}'
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
    # FIXME tweak legend if all the modifications are included
    plt.legend()
    mod_label = mod_name if mod_name else 'all'
    plt.savefig(f'{output_prefix}_{mod_label}.png')

def plot_PCA_all_points(pca, output_prefix):
    # add some annotations
    pca['mod_name'] = [s.split(' ')[0] for s in pca['sample name']]
    pca['concentration'] = [s.split(' ')[1] for s in pca['sample name']]
    print(output_prefix)
    # plot each of these, and all of them combined
    for mod in mod_names + [None]:
        print(mod)
        plot_PCA_for_mod(pca, output_prefix, mod)

def plot_all_sample_mean(PC_means, output_prefix):
    """Plots mean for all of the samples."""
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    for x_label, y_label in [('PC1','PC2')]:     # , ('PC2','PC3')]:
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

# plot_all_sample_mean(PC_mean_prism, 'PCA_prism')
# plot_all_sample_mean(PC_mean_scipy, 'PCA_scipy')

if True:
    for quantile_cutoff in [.01, .05, .1, .25, .5]:
        PC_scores = pandas.read_csv(f'../PCA_bases_5000reads/PCA_{quantile_cutoff}_quantile.csv.gz')
        PC_mean = PC_scores.groupby('sample name')[['PC1','PC2','PC3']].mean()
        PC_mean.reset_index(inplace=True)
        plot_all_sample_mean(PC_mean, f'{output_dir}/PCA_sample_means/PCA_quantile_cutoff_{quantile_cutoff}.png')

for quantile_cutoff in [.01, .05, .1, .25, .5]:
    PC_scores = pandas.read_csv(f'../PCA_bases_5000reads/PCA_{quantile_cutoff}_quantile.csv.gz')
    plot_PCA_all_points(PC_scores,
            f'{output_dir}/PCA_all_points/PCA_quantile_cutoff_{quantile_cutoff}')

