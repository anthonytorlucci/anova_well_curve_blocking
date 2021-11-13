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
import lasio
import binarytree

# local
from anova_log_blocking import anova_zoning, number_of_zones, _anova_recursive_tree_build
from log_curve_visualization import curve_plot

las_obj = lasio.read("data/us49025106110000_0_00028h489546.las")

df = las_obj.df()
df['DEPTH'] = df.index
df = df.fillna(method="ffill")
df = df.fillna(method="bfill")

depth = df.index.to_numpy()
curve = df['GRD'].to_numpy()

# select a smaller depth range
z_lower_indx = numpy.argwhere(depth > 725)
idx0 = z_lower_indx[0,0]
z_upper_indx = numpy.argwhere(depth < 875)
idx1 = z_upper_indx[-1,0]

depth = depth[idx0:idx1]
curve = curve[idx0:idx1]

root = _anova_recursive_tree_build(node=binarytree.Node(value=0), a=curve, min_samples_in_zone=25)
print(root)

# END