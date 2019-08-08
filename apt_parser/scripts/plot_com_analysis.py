#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.patches import Rectangle
import glob
import os
import apt_parser
import logging
apt_parser.init_log()

LOG = logging.getLogger('parse_apt')

folder = "../../xml"
files = glob.glob(os.path.join(folder, "*.aptx"))

(APT_sim, mem_volume, time_volume) = apt_parser.analyse_apt_list(files)

with open("commissioning_analysis.dat", 'w') as of:
    of.write("APT_sim = {\n")
    for k in sorted(list(APT_sim.keys())):
        v = APT_sim[k]
        of.write("'{}': {},\n".format(k,v))
    of.write("}\n")

    of.write("mem_volume = {\n")
    for k in sorted(list(mem_volume.keys())):
        v = mem_volume[k]

        formatted_list = ", ".join(map(lambda x: "{:.1f}".format(x), v))
        of.write("'{}': [{}]\n".format(k,formatted_list))
    of.write("}\n")
    of.write("time_volume = {\n")

    for k in sorted(list(time_volume.keys())):
        v = time_volume[k]
        formatted_list = ", ".join(map(lambda x: "{:.2f}".format(x), v))

        of.write("'{}': [{}]\n".format(k, formatted_list))
    of.write("}\n")

cars = list(time_volume.keys())
cars.sort()



total_time = np.zeros((len(cars)))
max_time   = np.zeros_like(total_time)
max_ram    = np.zeros_like(total_time)
total_time_dict   = {}
max_ram_dict    = {}
simulations =  []

for (i, car) in enumerate(cars):
    total_time[i] = time_volume[car].sum()
    max_time[i] = time_volume[car].sum()
    max_ram[i] = mem_volume[car].max()

    total_time_dict[car] = time_volume[car].sum()
    max_ram_dict[car] = mem_volume[car].max()

mem_threshold = 16.  # GB
# Split CARs in big and low memory
(big_mem_mask,) = np.where(max_ram > mem_threshold)


colors = iter(cm.nipy_spectral(np.linspace(0, 1, len(big_mem_mask))))



x = []
y = []
names = []
c = []

for i in big_mem_mask:
    car = cars[i]
    memory = mem_volume[car]  # Memory in GB
    time = time_volume[car] / 3600.  # Time in hours

    # We have to repeat a subarray for each exposure and dither
    name = []
    for d in APT_sim[car]:
        tmp = ["{}: {}, {}, {}".format(d["filename"], d["template"], d["obs_id"], d["subarray"])] * (d["exposures"] * d["NDither"])
        if d["subarray"] == "BOTH":
            tmp *= 2
        try:
            if d["detector"] == "ALL":
                tmp *= 3
        except KeyError:
            pass
        name.extend(tmp)
    color = next(colors)

    x.extend(memory)
    y.extend(time)
    c.extend([color]*len(time))
    names.extend(name)


def stacked_bar(bar_dict, ylabel=None, cl_label=None):
    """

    :param bar_dict: Each bar is a key in this dictionnary. The corresponding value is a numpy array containing all the Y-value for that bar
    :type bar_dict: dict(np.ndarray)
    :param ylabel: Unit of the Y axis
    :type ylabel: str
    :return:
    :rtype:
    """
    bar_w = 1
    bar_sep = 1
    edgecolor = "#000000"
    default_color = "#bbbbbb"
    linewidth = 0.5
    plot_margin = 1
    x_start = 0
    y_start = 0
    bar_labels = list(bar_dict.keys())
    bar_labels.sort()
    nbars = len(bar_labels)

    # We convert values to numpy array if necessary
    # flatten ensure that a single value is converted to a list.
    # Else we have a 0-d array error when trying to loop on it.
    for key in bar_labels:
        bar_dict[key] = np.asarray(bar_dict[key]).flatten()

    max_value = 0
    nb_values = 0
    for val in bar_dict.values():
        max_value = max(max_value, val.max())
        nb_values = max(nb_values, val.size)

    fig, ax = plt.subplots()
    fig.subplots_adjust(left=0.1, bottom=0.15, right=0.89, top=0.96, wspace=0.26, hspace=0.26)

    if nb_values > 1:
        cmap = mpl.cm.rainbow
        norm = mpl.colors.Normalize(vmin=0, vmax=max_value)

        # Get normalized color for given value
        get_color = lambda x: cmap(x/max_value)


        cbar_ax = fig.add_axes([0.9, 0.05, 0.025, 0.9])  # [left, bottom, width, height]
        cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap,
                                        norm=norm,
                                        orientation='vertical')
        if cl_label:
            cb1.set_label(cl_label)
    else:
        # We have only one value per bar, we define the same color everywhere
        get_color = lambda x: default_color

    x = x_start
    xmin = x - plot_margin
    ymin = y_start  # No margin for Ymin, to avoid ugly separation between 0 and the bottomline of the plot
    ymax = ymin
    tick_position = []

    for key in bar_labels:
        bar = bar_dict[key]
        height = y_start
        tick_position.append(x+bar_w/2)
        for val in bar:
            y = height
            height += val
            ymax = max(ymax, height) # We update ymax if necessary
            rect = Rectangle((x, y), bar_w , val , facecolor=get_color(val), edgecolor=edgecolor, linewidth=linewidth)
            ax.add_patch(rect)
        x += bar_w + bar_sep

    xmax = x - bar_sep + plot_margin
    ymax += plot_margin
    ax.set_xlim(xmin=xmin, xmax=xmax)
    ax.set_ylim(ymin=ymin, ymax=ymax)

    ax.set_xticks(tick_position)
    ax.set_xticklabels(bar_labels, rotation="vertical")

    if ylabel:
        ax.set_ylabel(ylabel)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, which='major', color='#000000', linestyle='--')


    return fig


def hover_plot(x, y, names, colors, xlabel, ylabel):
    """
    Scatter plot (X,Y) with description "names" added to a dot when hovering

    Use plt.show() to display the plot

    :param x: X axis values
    :type x: array(float)
    :param y: Y axis values
    :type y: array(float)
    :param names: hover description info for each dot
    :type names: list(str)
    :param colors: color for each dot
    :type colors: list
    :param xlabel: Label for X axis
    :type xlabel: str
    :param ylabel: label for Y-axis
    :type ylabel: str

    """


    norm = plt.Normalize(1,4)
    cmap = plt.cm.RdYlGn

    fig,ax = plt.subplots()
    sc = ax.scatter(x, y, c=colors, s=100, cmap=cmap, norm=norm)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(ind):

        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        # text = " ".join([names[n] for n in ind["ind"]])
        text = "\n".join([names[n] for n in ind["ind"]])
        # text = names[ind["ind"][0]]
        annot.set_text(text)
        annot.get_bbox_patch().set_facecolor(colors[ind["ind"][0]])
        annot.get_bbox_patch().set_alpha(0.4)


    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)

print("CAR & nb files & min RAM [GB] & Total time [h] \\\\\hline")
nb_exposures = 0
full_time = 0
for car in cars:
    msg = "{} & {} & {:.1f} GB & {:.1f} h\\\\\hline".format(car, len(mem_volume[car]), max_ram_dict[car], total_time_dict[car])
    print(msg)
    nb_exposures += len(mem_volume[car])
    full_time += total_time_dict[car]

print("Total number of exposures: {}".format(nb_exposures))
print("Total time for whole commissioning: {:.1f} h".format(full_time))

hover_plot(x, y, names, c, xlabel="Memory [GB]", ylabel="Pipeline Compute time [h]")
fig_time = stacked_bar(time_volume, ylabel="Total Compute time [h]", cl_label="Individual Compute time [h]")
fig_time.savefig("commissioning_time.pdf", bbox_inches="tight", pad_inches=0.1)

fig_max_ram = stacked_bar(max_ram_dict, ylabel="Max RAM used [GB]")
fig_max_ram.savefig("commissioning_max_memory.pdf", bbox_inches="tight", pad_inches=0.1)

plt.show()

# TODO plots barre cumulative pour chaque CAR (cumulé en temps, et en GB)
#TODO enumérer les problèmes du pipeline que mon paquet corrige (threads Numpy, effets de bords parce que le
# process n'est pas tué, gestion de la RAM)
#TODO faire 2e figure avec le fit sur l'autre ordinateur pour montrer la dépendance en puissance du CPU