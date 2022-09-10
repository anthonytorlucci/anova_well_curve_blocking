"""create a gif from images of anova log blocking with increasing minimum number of samples in zone"""

# standard libs
import os

# 3rd party
import numpy
import pandas
import matplotlib
import matplotlib.pyplot as plt
from PIL import Image

# local
from anova_log_blocking import anova_zoning, number_of_zones
from log_curve_visualization import curve_plot

df = pandas.read_csv('data/490252319500_GammaRay.csv')
df = df.set_index('DEPT')

depth = df.index.to_numpy()
curve = df['GR'].to_numpy()
print(f'number of zones in original curve: {number_of_zones(arr=curve)}')
min_depth = 2500
max_depth = 3000
z_lower_indx = numpy.argwhere(depth > min_depth)
idx0 = z_lower_indx[0,0]
z_upper_indx = numpy.argwhere(depth < max_depth)
idx1 = z_upper_indx[-1,0]
depth = depth[idx0:idx1]
curve = curve[idx0:idx1]
print(f'number of zones in depth range {min_depth} to {max_depth}: {number_of_zones(arr=curve)}')


nsamples = [n for n in range(2,48)]
for n in nsamples:
    curve_anova = anova_zoning(input_array=curve, min_samples_in_zone=n)
    # print(f'number of zones with {n} samples window: {number_of_zones(arr=curve_anova)}')
    # print(f'mean squared error with {n} samples window: {numpy.mean(numpy.square(curve - curve_anova))}')

    # create plot
    fig, axs = plt.subplots(ncols=2, nrows=1, constrained_layout=True)  # figsize=(16.5,9.5)
    fig.suptitle('constant thickness log blocking')
    curve_plot(ax=axs[0], depth=depth, curve=curve, ax_title='GR', cmap_name='gist_earth')
    axs[0].set_ylabel('DEPTH')
    curve_plot(ax=axs[1], depth=depth, curve=curve_anova, ax_title=f'n={n}', cmap_name='gist_earth')
    plt.savefig(fname=f'./images/anova_blocked_curve_{n:02}.png')

mpl_images = [f'./images/anova_blocked_curve_{n:02}.png' for n in range(2,48)]
frames = [Image.open(image) for image in mpl_images]
frame = frames[0]
frame.save("anova.gif", format="GIF", append_images=frames, save_all=True, duration=100, loop=0)