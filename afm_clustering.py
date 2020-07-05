#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# File      : afm_clustering.py
# Author    : Chris Jones
# Email     : crj341@student.bham.ac.uk
# Date      : 26/06/2020

# This script applies the HDBSCAN clustering algorithm to find clusters of
# data based on the adhesion vs. modulus plot. Data points outside any clusters
# are marked as noise. This method allowd for more accurate calculation of mean
# adhesion and modulus values, with low standard deviations. However, should be
# used with care to ensure you are not 'over cleaning' the data by removing
# valid data points. The hdbscan algorithm can be tuned by adjusting the
# min_cluster_size parameter, and others described in the liteature:
# https://hdbscan.readthedocs.io/en/latest/index.html
# The method may be useful for non uniform samples, allowing automated measurement
# of adhesion and modulus for different phases across the sample surface.

# Import AFM class
from    afm import  AFM

# Import packages for clustering
import  hdbscan
import  numpy               as  np
import  matplotlib.pyplot   as  plt
import  seaborn             as  sns
import  sklearn.datasets    as  data

# Run AFM functions
expt_1 = AFM()
expt_1.run()

# Shape adhesion and modulus data
_adhesion = expt_1.adhesion
_modulus = expt_1.modulus
_adhesion = _adhesion.reshape(-1,1)
_modulus = _modulus.reshape(-1,1)
data = np.append(_adhesion, _modulus, axis = 1)

# Find clusters
clusterer = hdbscan.HDBSCAN(
    allow_single_cluster = True,
    min_cluster_size = 10
)
clusterer.fit(data)

# Plot condensed tree
clusterer.condensed_tree_.plot(
    select_clusters=True,
    selection_palette=sns.color_palette()
)

# Plot clusters
palette = sns.color_palette()
cluster_colors = [sns.desaturate(palette[col], sat)
                  if col >= 0 else (0.5, 0.5, 0.5) for col, sat in
                  zip(clusterer.labels_, clusterer.probabilities_)]

plt.scatter(_adhesion, _modulus, c=cluster_colors)

# Analyse results
mean_adhesion = []
std_adhesion = []
mean_modulus = []
std_modulus = []

for i in np.unique(clusterer.labels_):
    if i < 0:
        continue

    mean_adhesion = np.append(mean_adhesion, np.mean(_adhesion[clusterer.labels_ == i]))
    std_adhesion = np.append(std_adhesion, np.std(_adhesion[clusterer.labels_ == i]))
    mean_modulus = np.append(mean_modulus, np.mean(_modulus[clusterer.labels_ == i]))
    std_modulus = np.append(std_modulus, np.std(_modulus[clusterer.labels_ == i]))

print('For each cluster:')
print('Mean adhesion:', mean_adhesion)
print('St. dev.  adhesion:', std_adhesion)
print('Mean modulus:', mean_modulus)
print('St. dev.  modulus:', std_modulus)