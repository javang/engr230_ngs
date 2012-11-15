
import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib as mpl
from math import pi


class MyColors:
    colors_list = ['red', 'blue', 'green', 'orange', 'purple', 'gray',
               'orange', 'DarkOrchid', 'olive', 'salmon','MediumPurple','LightSteelBlue',
               'LightCoral', 'gold', 'DarkSlateGray', 'fuchsia', 'forestgreen', 'yellow','pink']
    def __init__(self):
        self.i = 0

    def get_next_color(self):
        self.i += 1
        if self.i == len(self.colors_list):
            self.i = 0
        return self.colors_list[self.i]

def fig2(coverages, gc_contents, lengths, genuses, fn_output):

    ############# SETUP FIGURE ####################
    plt.rcParams["axes.linewidth"] = 1
    dpi = 900
    fontsize = 6
    legend_fontsize = 5
    fig_height = 4
    fig_width = 8
    figure = plt.figure(figsize=(fig_width, fig_height), dpi=dpi, facecolor='w', edgecolor='k')
    # axes
    ax = figure.add_subplot(111)
    figure.sca(ax)
    ax.grid()
    ax.set_xlabel("Coverage", fontsize=fontsize)
    ax.set_xscale("log")
    ax.set_xlim(left=1,right=max(coverages))
    ax.set_ylabel("% CG", fontsize=fontsize)
    ax.set_position([0.05,0.1,0.65,0.8])
    plt.setp(ax.get_yticklabels(), fontsize=fontsize)
    plt.setp(ax.get_xticklabels(), fontsize=fontsize)
    # set the scale so the longest scaffold uses 1/100000 of the pixels in the figure
    total_pixels = dpi**2 * fig_height * fig_width
    scatter_scale = 1e-5 * total_pixels / max(lengths)

    # if not assigned
    if not genuses:
        sizes = np.array(lengths) * scatter_scale
        sc = ax.scatter(coverages, gc_contents, s=sizes, c='gray',
                        marker='o', lw=None, edgecolor="none")
        plt.savefig(fn_output,dpi=dpi)
        return

    # dictionary of indices with the name of the genus as key
    unique_genuses = set(genuses)
    nomatch = "No match"
    mycolors = MyColors()
    scatter_plots = []
    labels = []
    for ug in unique_genuses:
        if ug == nomatch:
            color = 'gray'
        else:
            color = mycolors.get_next_color()
        coverages_to_plot = []
        gc_contents_to_plot = []
        lengths_to_plot = []
        for i, g in enumerate(genuses):
            if ug == g:
                coverages_to_plot.append(coverages[i])
                coverages.pop(i)
                gc_contents_to_plot.append(gc_contents[i])
                gc_contents.pop(i)
                genuses.pop(i)
                lengths_to_plot.append(lengths[i])
                lengths.pop(i)
        sizes = np.array(lengths_to_plot) * scatter_scale
        sc = ax.scatter(coverages_to_plot, gc_contents_to_plot, s=sizes, c=color,
                        marker='o', lw=None, edgecolor="none")
        labels.append(ug)
        scatter_plots.append(sc)

    plt.figlegend(scatter_plots,labels,
                  #bbox_to_anchor=(1.05, 1), # put slightly to the rigth of the plot
                  loc="upper right",
                  ncol=2,
                  prop={"size":legend_fontsize})
    plt.savefig(fn_output,dpi=dpi)

