"""log_curve_visualization.py
utility functions to plot log curve data
"""

# standard libs

# 3rd party
import matplotlib
import matplotlib.pyplot as plt
import numpy

# local

# "Enhancing Visualization of Well Logs With Plot Fills"
# https://towardsdatascience.com/enhancing-visualization-of-well-logs-with-plot-fills-72d9dcd10c1b
def curve_plot(ax:matplotlib.axes.Axes, depth:numpy.ndarray, curve:numpy.ndarray, ax_title='', cmap_name='nipy_spectral'):
    # setting up the viz
    left_col_value = 0
    right_col_value = 150

    #calculate the span of values
    span = abs(left_col_value - right_col_value)

    #assign a color map
    cmap = plt.get_cmap(cmap_name)

    #create array of values to divide up the area under curve
    color_index = numpy.arange(left_col_value, right_col_value, span / 100)

    ax.set_ylim(depth[-1], depth[0])
    # ax.set_ylim(depth[-1], depth[-101])
    ax.set_xlim(left_col_value, right_col_value)
    ax.set_title(ax_title)
    ax.grid(True)

    #loop through each value in the color_index
    for index in sorted(color_index):
        index_value = (index - left_col_value)/span
        color = cmap(index_value) #obtain colour for color index value
        ax.fill_betweenx(depth, 0 , curve, where = curve >= index,  color = color)