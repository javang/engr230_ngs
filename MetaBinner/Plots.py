import matplotlib
# set a backend that does not require interactive mode. Useful for remote figure generation
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#if plt.isinteractive():
#    plt.ioff()
import matplotlib.cm as cm
import matplotlib as mpl
import numpy as np
from math import pi

import MetaBinner.definitions as defs

import csv
import sys


"""
# dictionary of colors to make the figures uniform. For some reason matplotlib
# does not respect the order
Genus2Color ={ "thermobaculum": "red",
               "thermomicrobium": "blue",
               "sphaerobacter": "yellow",
                "patulibacter": "indigo",
                "nonomuraea": "salmon",
                "conexibacter": "cyan",
                "thermobacillus": "gold",
                "paenibacillus": "orange",
                "cohnella": "khaki",
                "caldalkalibacillus": "lime",
                "sulfurihydrogenibium": "lightcoral",
                "rhodothermus": "darkorchid",
                "hyphomicrobium": "green",
                "gemmatimonas": "mediumpurple",
                "thermus": "pink",
                "thermobispora": "black"
                }
"""

#,'LightSteelBlue'
class MyColors:
    colors_list = ['red', 'blue', 'green', 'orange', 'yellow', 'purple',
                'DarkOrchid', 'olive', 'salmon','brown','MediumPurple',
               'LightCoral', 'gold', 'pink','cyan', 'BurlyWood',
               'chocolate', 'crimson','deeppink', 'black','aquamarine','aqua',
               'Darkorange','darkmagenta','lightgreen', 'greenyellow', "indigo", "lime",
               "MidnightBlue", "SaddleBrown",'indianred','khaki']
    def __init__(self):
        self.i = 0

    def get_next_color(self):
        self.i += 1
        if self.i == len(self.colors_list):
            self.i = 0
        return self.colors_list[self.i]

def fig2(coverages, gc_contents, lengths, genera, fn_output):

    ############# SETUP FIGURE ####################
    plt.rcParams["axes.linewidth"] = 1
    dpi = 900
    fontsize = 6
    legend_fontsize = 5
    fig_height = 4
    fig_width = 8
    figure = plt.figure(figsize=(fig_width, fig_height),
                                    dpi=dpi, facecolor='w', edgecolor='k')
    # axes
    ax = figure.add_subplot(111)
    figure.sca(ax)
    ax.grid()

    ax.set_xlabel("Coverage", fontsize=fontsize)
    ax.set_xscale("log")
    ax.set_xlim(left=1,right=max(coverages))
    ax.set_ylabel("CG fraction", fontsize=fontsize)

#    ax.set_xlabel("PCA1", fontsize=fontsize) # change the labels for PCA Plots
#    ax.set_ylabel("PCA2", fontsize=fontsize)

    ax.set_position([0.07,0.1,0.65,0.8])
    plt.setp(ax.get_yticklabels(), fontsize=fontsize)
    plt.setp(ax.get_xticklabels(), fontsize=fontsize)
    # set the scale so the longest scaffold uses 1/100000 of the pixels in the figure
    total_pixels = dpi**2 * fig_height * fig_width
    scatter_scale = 8e-6 * total_pixels / max(lengths)

    # if there is no data for the assignments, plot the raw figure
    if not genera:
        sizes = np.array(lengths) * scatter_scale
        sc = ax.scatter(coverages, gc_contents, s=sizes, c='gray',
                        marker='o', lw=None, edgecolor="none")
        plt.savefig(fn_output,dpi=dpi)
        return

    unique_genera = set(genera)
    if defs.not_assigned in unique_genera:
        unique_genera.remove(defs.not_assigned)

    unique_genera = [g for g in unique_genera]
    unique_genera.sort()
    unique_genera = [defs.not_assigned] + unique_genera
    mycolors = MyColors()
    scatter_plots = []
    labels = []
    for ug in unique_genera:
        if ug == defs.not_assigned:
            color = 'gainsboro'
        else:
            color = mycolors.get_next_color() # Genus2Color[ug]
        coverages_to_plot = []
        gc_contents_to_plot = []
        lengths_to_plot = []
        indices_to_pop = []
        for i, g in enumerate(genera):
            if ug == g:
                coverages_to_plot.append(coverages[i])
                gc_contents_to_plot.append(gc_contents[i])
                lengths_to_plot.append(lengths[i])
                indices_to_pop.append(i)

        # pop after finishing adding elements to avoid messing with the indices:
        indices_to_pop.sort()
        indices_to_pop.reverse()
        for i in indices_to_pop:
            coverages.pop(i)
            gc_contents.pop(i)
            genera.pop(i)
            lengths.pop(i)
        # if there is nothing to plot, continue
        if len(coverages_to_plot) == 0:
            continue
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


def fig3(coverages, gc_contents, lengths, genera, fn_output):
    # plot the graph of covgerage vs (coverage/GC content)
    ############# SETUP FIGURE ####################
    plt.rcParams["axes.linewidth"] = 1
    dpi = 900
    fontsize = 6
    legend_fontsize = 5
    fig_height = 4
    fig_width = 8

    figure = plt.figure(figsize=(fig_width, fig_height),
                                    dpi=dpi, facecolor='w', edgecolor='k')
    # axes
    ax = figure.add_subplot(111)
    figure.sca(ax)
    ax.grid()

    ax.set_xlabel("Coverage", fontsize=fontsize)
#    ax.set_xscale("log")
    ax.set_xlim(left=1,right=200)
    ax.set_ylabel("CG fraction", fontsize=fontsize)
#    ax.set_yscale("log")

#    ax.set_xlabel("PCA1", fontsize=fontsize) # change the labels for PCA Plots
#    ax.set_ylabel("PCA2", fontsize=fontsize)

    ax.set_position([0.07,0.1,0.65,0.8])
    plt.setp(ax.get_yticklabels(), fontsize=fontsize)
    plt.setp(ax.get_xticklabels(), fontsize=fontsize)
    # set the scale so the longest scaffold uses 1/100000 of the pixels in the figure
    total_pixels = dpi**2 * fig_height * fig_width
    scatter_scale = 8e-6 * total_pixels / max(lengths)

    # if there is no data for the assignments, plot the raw figure
    if not genera:
        sizes = np.array(lengths) * scatter_scale
        sc = ax.scatter(coverages, gc_contents, s=sizes, c='gray',
                        marker='o', lw=None, edgecolor="none")
        plt.savefig(fn_output,dpi=dpi)
        return

    unique_genera = set(genera)
    if defs.not_assigned in unique_genera:
        unique_genera.remove(defs.not_assigned)

    unique_genera = [g for g in unique_genera]
    unique_genera.sort()
    unique_genera = [defs.not_assigned] + unique_genera
    mycolors = MyColors()
    scatter_plots = []
    labels = []
    for ug in unique_genera:
        if ug == defs.not_assigned:
            color = 'gainsboro'
        else:
            color = mycolors.get_next_color() # Genus2Color[ug]
        coverages_to_plot = []
        gc_contents_to_plot = []
        lengths_to_plot = []
        indices_to_pop = []
        for i, g in enumerate(genera):
            if ug == g:
                coverages_to_plot.append(coverages[i])
                gc_contents_to_plot.append(gc_contents[i])
                lengths_to_plot.append(lengths[i])
                indices_to_pop.append(i)

        # pop after finishing adding elements to avoid messing with the indices:
        indices_to_pop.sort()
        indices_to_pop.reverse()
        for i in indices_to_pop:
            coverages.pop(i)
            gc_contents.pop(i)
            genera.pop(i)
            lengths.pop(i)
        # if there is nothing to plot, continue
        if len(coverages_to_plot) == 0:
            continue
        sizes = np.array(lengths_to_plot) * scatter_scale

        c = np.array(coverages_to_plot) 
        g = np.array(gc_contents_to_plot) 
        sc = ax.scatter(c, g, s=sizes, c=color,
                        marker='o', lw=None, edgecolor="none")
        labels.append(ug)
        scatter_plots.append(sc)

    plt.figlegend(scatter_plots,labels,
                  #bbox_to_anchor=(1.05, 1), # put slightly to the rigth of the plot
                  loc="upper right",
                  ncol=2,
                  prop={"size":legend_fontsize})
    plt.savefig(fn_output,dpi=dpi)

