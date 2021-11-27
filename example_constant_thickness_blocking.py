"""example_constant_thickness_blocking.py"""

# standard libs
import os

# 3rd party
import numpy
import pandas
import matplotlib
import matplotlib.pyplot as plt

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

df = pandas.read_csv('data/490252319500_GammaRay.csv')
df = df.set_index('DEPT')

depth = df.index.to_numpy()
curve = df['GR'].to_numpy()

samples = [3,6,12,24]
blocked_curves = []

print("number of zones in original curve: ", number_of_zones(arr=curve))
for i in range(len(samples)):
    curve_thick = constant_thickness_zoning(input_array=curve, nsamples=samples[i])
    print("number of zones with {} samples window: {}".format(samples[i], number_of_zones(arr=curve_thick)))
    print("mean squared error with {} samples window: {}".format(samples[i], numpy.mean(numpy.square(curve - curve_thick))))
    blocked_curves.append(curve_thick)

# create plot
gs_kw = dict(width_ratios=[1,1,1,1,1], height_ratios=[2,1])
fig, axs = plt.subplots(ncols=len(samples)+1, nrows=2, constrained_layout=True, gridspec_kw=gs_kw, figsize=(16.5,9.5))
fig.suptitle('constant thickness log blocking')
curve_plot(ax=axs[0,0], depth=depth, curve=curve, ax_title='GR', cmap_name='gist_earth')
axs[0,0].set_ylabel('DEPTH')
for i in range(len(samples)):
    c = blocked_curves[i]
    curve_plot(ax=axs[0,i+1], depth=depth, curve=c, ax_title='nsamp={}'.format(samples[i]), cmap_name='gist_earth')

# plot depth range
z_lower_indx = numpy.argwhere(depth > 2750)
idx0 = z_lower_indx[0,0]
z_upper_indx = numpy.argwhere(depth < 2950)
idx1 = z_upper_indx[-1,0]
curve_plot(ax=axs[1,0], depth=depth[idx0:idx1], curve=curve[idx0:idx1], ax_title='GR', cmap_name='gist_earth')
axs[1,0].set_ylabel('DEPTH')
for i in range(len(samples)):
    c = blocked_curves[i]
    curve_plot(ax=axs[1,i+1], depth=depth[idx0:idx1], curve=c[idx0:idx1], ax_title='nsamp={}'.format(samples[i]), cmap_name='gist_earth')

#plt.show()
plt.savefig('images/constant_thickness_blocking.png')