"""example_tree.py
prints the binary tree structure
"""

# standard libs
import os

# 3rd party
import numpy
import pandas
import matplotlib
import matplotlib.pyplot as plt
import binarytree

# local
from anova_log_blocking import anova_zoning, number_of_zones, _anova_recursive_tree_build
from log_curve_visualization import curve_plot

df = pandas.read_csv('data/490252319500_GammaRay.csv')
df = df.set_index('DEPT')
#++print(df)

# if missing values, use forward fill and backward fill
if df.isnull().values.any():
    df = df.fillna(method="ffill")
    df = df.fillna(method="bfill")

depth = df.index.to_numpy()
curve = df['GR'].to_numpy()

# select a smaller depth range
z_lower_indx = numpy.argwhere(depth > 2800)
idx0 = z_lower_indx[0,0]
z_upper_indx = numpy.argwhere(depth < 2900)
idx1 = z_upper_indx[-1,0]

depth = depth[idx0:idx1]
curve = curve[idx0:idx1]

root = _anova_recursive_tree_build(node=binarytree.Node(value=0), a=curve, min_samples_in_zone=6)
print(root)

# END