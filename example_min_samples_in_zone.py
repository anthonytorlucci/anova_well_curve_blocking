"""example_min_samples_in_zone.py"""

# standard libs
import os

# 3rd party
import numpy
import pandas
import matplotlib
import matplotlib.pyplot as plt
import lasio

# local
from anova_log_blocking import anova_zoning, number_of_zones
from log_curve_visualization import curve_plot

las_obj = lasio.read("data/us49025106110000_0_00028h489546.las")
# print(las_obj.well)
# print(las_obj.curves)
# print(type(las_obj))

df = las_obj.df()
df['DEPTH'] = df.index
df = df.fillna(method="ffill")
df = df.fillna(method="bfill")

depth = df.index.to_numpy()
curve = df['GRD'].to_numpy()

print("number of zones in original curve: ", number_of_zones(arr=curve))
curve_anova05 = anova_zoning(input_array=curve, min_samples_in_zone=5)
print("number of zones with min 5 samples window: ", number_of_zones(arr=curve_anova05))
print("mean squared error with min 5 samples window: ", numpy.mean(numpy.square(curve - curve_anova05)))
curve_anova10 = anova_zoning(input_array=curve, min_samples_in_zone=10)
print("number of zones with min 10 samples window: ", number_of_zones(arr=curve_anova10))
print("mean squared error with min 10 samples window: ", numpy.mean(numpy.square(curve - curve_anova10)))
curve_anova15 = anova_zoning(input_array=curve, min_samples_in_zone=15)
print("number of zones with min 15 samples window: ", number_of_zones(arr=curve_anova15))
print("mean squared error with min 15 samples window: ", numpy.mean(numpy.square(curve - curve_anova15)))
curve_anova20 = anova_zoning(input_array=curve, min_samples_in_zone=20)
print("number of zones with min 20 samples window: ", number_of_zones(arr=curve_anova20))
print("mean squared error with min 20 samples window: ", numpy.mean(numpy.square(curve - curve_anova20)))

# create plot
gs_kw = dict(width_ratios=[1,1,1,1,1], height_ratios=[2,1])
fig, axs = plt.subplots(ncols=5, nrows=2, constrained_layout=True, gridspec_kw=gs_kw, figsize=(16.5,9.5))
fig.suptitle('anova log blocking')
curve_plot(ax=axs[0,0], depth=depth, curve=curve, ax_title='GRD', cmap_name='gist_earth')
axs[0,0].set_ylabel('DEPTH')
curve_plot(ax=axs[0,1], depth=depth, curve=curve_anova05, ax_title='min_samp=5', cmap_name='gist_earth')
curve_plot(ax=axs[0,2], depth=depth, curve=curve_anova10, ax_title='min_samp=10', cmap_name='gist_earth')
curve_plot(ax=axs[0,3], depth=depth, curve=curve_anova15, ax_title='min_samp=15', cmap_name='gist_earth')
curve_plot(ax=axs[0,4], depth=depth, curve=curve_anova20, ax_title='min_samp=20', cmap_name='gist_earth')
# plot depth range
z_lower_indx = numpy.argwhere(depth > 725)
# print(z_lower_indx.shape)
# print(z_lower_indx[0,0])
idx0 = z_lower_indx[0,0]
z_upper_indx = numpy.argwhere(depth < 875)
# print(z_upper_indx.shape)
# print(z_upper_indx[-1,0])
idx1 = z_upper_indx[-1,0]
curve_plot(ax=axs[1,0], depth=depth[idx0:idx1], curve=curve[idx0:idx1], ax_title='GRD', cmap_name='gist_earth')
axs[1,0].set_ylabel('DEPTH')
curve_plot(ax=axs[1,1], depth=depth[idx0:idx1], curve=curve_anova05[idx0:idx1], ax_title='min_samp=5', cmap_name='gist_earth')
curve_plot(ax=axs[1,2], depth=depth[idx0:idx1], curve=curve_anova10[idx0:idx1], ax_title='min_samp=10', cmap_name='gist_earth')
curve_plot(ax=axs[1,3], depth=depth[idx0:idx1], curve=curve_anova15[idx0:idx1], ax_title='min_samp=15', cmap_name='gist_earth')
curve_plot(ax=axs[1,4], depth=depth[idx0:idx1], curve=curve_anova20[idx0:idx1], ax_title='min_samp=20', cmap_name='gist_earth')
# plt.show()
plt.savefig('images/anova_min_samples_in_zone.png')