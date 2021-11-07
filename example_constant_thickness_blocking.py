"""example_constant_thickness_blocking.py"""

# standard libs
import os

# 3rd party
import numpy
import pandas
import matplotlib
import matplotlib.pyplot as plt
import lasio

# local
from anova_log_blocking import number_of_zones
from log_curve_visualization import curve_plot

def constant_thickness_zoning(input_array:numpy.ndarray, nsamples:int):
    output_array = numpy.zeros_like(input_array)
    
    for n in range(0,len(input_array),nsamples):
        output_array[n:n+nsamples] = numpy.median(input_array[n:n+nsamples])
    # handle the last samples
    if n+nsamples < len(input_array):
        output_array[n+nsamples:] = numpy.median(input_array[n+nsamples:])
    return output_array

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
curve_thick05 = constant_thickness_zoning(input_array=curve, nsamples=5)
print("number of zones with 5 samples window: ", number_of_zones(arr=curve_thick05))
print("mean squared error with 5 samples window: ", numpy.mean(numpy.square(curve - curve_thick05)))
curve_thick10 = constant_thickness_zoning(input_array=curve, nsamples=10)
print("number of zones with 10 samples window: ", number_of_zones(arr=curve_thick10))
print("mean squared error with 10 samples window: ", numpy.mean(numpy.square(curve - curve_thick10)))
curve_thick15 = constant_thickness_zoning(input_array=curve, nsamples=15)
print("number of zones with 15 samples window: ", number_of_zones(arr=curve_thick15))
print("mean squared error with 15 samples window: ", numpy.mean(numpy.square(curve - curve_thick15)))
curve_thick20 = constant_thickness_zoning(input_array=curve, nsamples=20)
print("number of zones with 20 samples window: ", number_of_zones(arr=curve_thick20))
print("mean squared error with 20 samples window: ", numpy.mean(numpy.square(curve - curve_thick20)))

# create plot
gs_kw = dict(width_ratios=[1,1,1,1,1], height_ratios=[2,1])
fig, axs = plt.subplots(ncols=5, nrows=2, constrained_layout=True, gridspec_kw=gs_kw, figsize=(16.5,9.5))
fig.suptitle('constant thickness log blocking')
curve_plot(ax=axs[0,0], depth=depth, curve=curve, ax_title='GRD', cmap_name='gist_earth')
axs[0,0].set_ylabel('DEPTH')
curve_plot(ax=axs[0,1], depth=depth, curve=curve_thick05, ax_title='nsamp=5', cmap_name='gist_earth')
curve_plot(ax=axs[0,2], depth=depth, curve=curve_thick10, ax_title='nsamp=10', cmap_name='gist_earth')
curve_plot(ax=axs[0,3], depth=depth, curve=curve_thick15, ax_title='nsamp=15', cmap_name='gist_earth')
curve_plot(ax=axs[0,4], depth=depth, curve=curve_thick20, ax_title='nsamp=20', cmap_name='gist_earth')
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
curve_plot(ax=axs[1,1], depth=depth[idx0:idx1], curve=curve_thick05[idx0:idx1], ax_title='nsamp=5', cmap_name='gist_earth')
curve_plot(ax=axs[1,2], depth=depth[idx0:idx1], curve=curve_thick10[idx0:idx1], ax_title='nsamp=10', cmap_name='gist_earth')
curve_plot(ax=axs[1,3], depth=depth[idx0:idx1], curve=curve_thick15[idx0:idx1], ax_title='nsamp=15', cmap_name='gist_earth')
curve_plot(ax=axs[1,4], depth=depth[idx0:idx1], curve=curve_thick20[idx0:idx1], ax_title='nsamp=20', cmap_name='gist_earth')
# plt.show()
plt.savefig('images/constant_thickness_blocking.png')