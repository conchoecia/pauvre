#!/usr/bin/env python
# pauvre - just a pore PhD student's plotting package
# Copyright (c) 2016-2017 Darrin T. Schultz. All rights reserved.
#
# This file is part of pauvre.
#
# pauvre is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pauvre is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pauvre.  If not, see <http://www.gnu.org/licenses/>.

import matplotlib
import platform
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import os.path as opath
from pauvre.functions import parse_fastq_length_meanqual
from pauvre.stats import stats
import pauvre.rcparams as rc
import logging

if platform.system() == 'Linux':
    matplotlib.use('agg')


# logging
logger = logging.getLogger('pauvre')


def generate_length_histogram_panel():
    return


def margin_plot(args):
    rc.update_rcParams()
    length, mean_qual = parse_fastq_length_meanqual(args.fastq)
    stats(args.fastq, length, mean_qual)

    if args.maxlen:
        max_plot_length = args.maxlen
    else:
        max_plot_length = int(np.percentile(length, 99))

    if args.lengthbin:
        length_bin_interval = args.lengthbin
    else:
        # Dividing by 80 is based on what looks good from experience
        length_bin_interval = int(max_plot_length / 80)

    length_bins = np.arange(0, max_plot_length, length_bin_interval)

    if args.maxqual:
        max_plot_qual = args.maxqual
    else:
        max_plot_qual = max(np.ceil(mean_qual))
    if args.qualbin:
        qual_bin_interval = args.qualbin
    else:
        # again, this is just based on what looks good from experience
        qual_bin_interval = max_plot_qual / 85
    qual_bins = np.arange(0, max_plot_qual, qual_bin_interval)

    # 250, 231, 34 light yellow
    # 67, 1, 85
    # R=np.linspace(65/255,1,101)
    # G=np.linspace(0/255, 231/255, 101)
    # B=np.linspace(85/255, 34/255, 101)

    # R=65/255, G=0/255, B=85/255
    Rf = 65 / 255
    Bf = 85 / 255
    pdict = {'red': ((0.0, Rf, Rf),
                     (1.0, Rf, Rf)),
             'green': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
             'blue': ((0.0, Bf, Bf),
                      (1.0, Bf, Bf)),
             'alpha': ((0.0, 0.0, 0.0),
                       (1.0, 1.0, 1.0))
             }
    # Now we will use this example to illustrate 3 ways of
    # handling custom colormaps.
    # First, the most direct and explicit:
    purple1 = LinearSegmentedColormap('Purple1', pdict)

    # make the pandas dataset to query
    df = pd.DataFrame(list(zip(length, mean_qual)), columns=['length', 'meanQual'])
    # only keep the dataframes that are finite
    df = df.dropna()

    # set the figure dimensions
    fig_width = 1.61 * 3
    fig_height = 1 * 3
    fig = plt.figure(figsize=(fig_width, fig_height))

    # set the panel dimensions
    heat_map_panel_width = fig_width * 0.5
    heat_map_panel_height = heat_map_panel_width * 0.62

    # find the margins to center the panel in figure
    left_margin = bottom_margin = (1 / 6)

    # lengthPanel
    length_panel_width = (1 / 8)

    # the color Bar parameters
    legend_panel_width = (1 / 24)

    # define padding
    h_padding = 0.02
    v_padding = 0.05

    panels = []

    # Quality histogram panel
    qual_panel_left = left_margin + length_panel_width + h_padding
    qual_panel_bottom = bottom_margin
    qual_panel_width = heat_map_panel_width / fig_width
    qual_panel_height = length_panel_width * fig_width / fig_height
    qual_panel_rec = [qual_panel_left, qual_panel_bottom,
                      qual_panel_width, qual_panel_height]

    qual_panel = plt.axes(qual_panel_rec)

    qual_panel.tick_params(axis='both', which='both',
                           bottom='on', labelbottom='on',
                           left='on', labelleft='on',
                           right='off', labelright='off',
                           top='off', labeltop='off')

    panels.append(qual_panel)

    # Length histogram panel
    length_panel_left = left_margin
    length_panel_bottom = bottom_margin + qual_panel_height + v_padding
    length_panel_height = heat_map_panel_height / fig_height
    length_panel_rec = [length_panel_left, length_panel_bottom,
                        length_panel_width, length_panel_height]

    length_panel = plt.axes(length_panel_rec)

    length_panel.tick_params(axis='both', which='both',
                             bottom='on', labelbottom='on',
                             left='on', labelleft='on',
                             right='off', labelright='off',
                             top='off', labeltop='off')
    panels.append(length_panel)

    # Heat map panel
    heat_map_panel_left = left_margin + length_panel_width + h_padding
    heat_map_panel_bottom = bottom_margin + qual_panel_height + v_padding
    heat_map_panel_rec = [heat_map_panel_left, heat_map_panel_bottom,
                          heat_map_panel_width / fig_width,
                          heat_map_panel_height / fig_height]
    heat_map_panel = plt.axes(heat_map_panel_rec)

    heat_map_panel.tick_params(axis='both', which='both',
                               bottom='off', labelbottom='off',
                               left='off', labelleft='off',
                               right='off', labelright='off',
                               top='off', labeltop='off')

    panels.append(heat_map_panel)
    heat_map_panel.set_title(args.title)

    # Legend panel
    legend_panel_left = left_margin + length_panel_width + \
                        heat_map_panel_width / fig_width + h_padding
    legend_panel_bottom = bottom_margin + qual_panel_height + v_padding
    legend_panel_height = heat_map_panel_height / fig_height
    legend_panel_rec = [legend_panel_left, legend_panel_bottom,
                        legend_panel_width, legend_panel_height]
    legend_panel = plt.axes(legend_panel_rec)

    legend_panel.tick_params(axis='both', which='both',
                             bottom='off', labelbottom='off',
                             left='off', labelleft='off',
                             right='on', labelright='on',
                             top='off', labeltop='off')
    panels.append(legend_panel)

    # plot the length histogram on y-axis
    length_bins_values, bins2 = np.histogram(length, length_bins)
    length_panel.set_ylim([0, max_plot_length])
    length_panel.set_xlim([0, max(length_bins_values * 1.1)])

    for step in np.arange(0, len(length_bins_values), 1):
        left = 1
        bottom = length_bins[step]
        width = length_bins_values[step]
        height = length_bins[step + 1] - length_bins[step]

        rectangle1 = mplpatches.Rectangle((left, bottom), width, height,
                                          linewidth=0.0,
                                          facecolor=(0.5, 0.5, 0.5),
                                          edgecolor=(0, 0, 0))
        length_panel.add_patch(rectangle1)
    length_panel.spines['top'].set_visible(False)
    length_panel.spines['right'].set_visible(False)
    length_panel.spines['bottom'].set_visible(True)
    length_panel.set_ylabel('Read Length')

    # plot the length histogram on x-axis
    qual_bins_values, bins2 = np.histogram(mean_qual, qual_bins)
    qual_panel.set_ylim([0, max(qual_bins_values * 1.1)])
    qual_panel.set_xlim([0, max_plot_qual])
    for step in np.arange(0, len(qual_bins_values), 1):
        left = qual_bins[step]
        bottom = 0
        width = qual_bins[step + 1] - qual_bins[step]
        height = qual_bins_values[step]

        rectangle1 = mplpatches.Rectangle((left, bottom), width, height,
                                          linewidth=0.0,
                                          facecolor=(0.5, 0.5, 0.5),
                                          edgecolor=(0, 0, 0))
        qual_panel.add_patch(rectangle1)
    qual_panel.spines['top'].set_visible(False)
    qual_panel.spines['right'].set_visible(False)
    qual_panel.spines['left'].set_visible(True)
    qual_panel.set_xlabel('Phred Quality')
    qual_panel.set_ylabel('Count')

    # plot the length histogram on x-axis
    hex_this = df.query('length<{} and meanQual<{}'.format(
        max_plot_length, max_plot_qual))
    # print(hexThis)

    heat_map_panel.set_xlim([0, max_plot_qual])
    heat_map_panel.set_ylim([0, max_plot_length])
    # This single line controls plotting the hex bins in the panel
    hex_vals = heat_map_panel.hexbin(hex_this['meanQual'], hex_this['length'], gridsize=49,
                            linewidths=0.0, cmap=purple1)
    counts = hex_vals.get_array()

    # plot the colorbar
    # completely custom for more control
    legend_panel.set_xlim([0, 1])
    legend_panel.set_ylim([0, 1000])
    legend_panel.set_yticks([int(x) for x in np.linspace(0, 1000, 6)])
    legend_panel.set_yticklabels([int(x) for x in np.linspace(0, max(counts), 6)])
    for i in np.arange(0, 1001, 1):
        rgba = purple1(i / 1001)
        alpha = rgba[-1]
        facec = rgba[0:3]
        rectangle1 = mplpatches.Rectangle((0, i), 1, 1,
                                          linewidth=0.0,
                                          facecolor=facec,
                                          edgecolor=(0, 0, 0),
                                          alpha=alpha)
        legend_panel.add_patch(rectangle1)
    legend_panel.spines['top'].set_visible(False)
    legend_panel.spines['left'].set_visible(False)
    legend_panel.spines['bottom'].set_visible(False)
    legend_panel.yaxis.set_label_position("right")
    legend_panel.set_ylabel('Number of Reads')

    for each in heat_map_panel.spines:
        heat_map_panel.spines[each].set_visible(False)

    file_base = opath.splitext(opath.basename(args.fastq))[0]

    for each in args.fileform:
        out_name = "{}.{}".format(file_base, each)
        if each == 'png':

            # transparent parameter not functioning. Need to fix.
            plt.savefig(out_name, dpi=args.dpi, transparent=False)
        else:
            plt.savefig(out_name, transparent=args.transparent)


def run(args):
    margin_plot(args)
    print(args.transparent)
