#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

import ast
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import os.path as opath
from sys import stderr
from pauvre.functions import print_images
from pauvre.stats import stats
import pauvre.rcparams as rc
import sys
import logging

# logging
logger = logging.getLogger('pauvre')


def generate_panel(panel_left, panel_bottom, panel_width, panel_height,
                   axis_tick_param='both', which_tick_param='both',
                   bottom_tick_param='on', label_bottom_tick_param='on',
                   left_tick_param='on', label_left_tick_param='on',
                   right_tick_param='off', label_right_tick_param='off',
                   top_tick_param='off', label_top_tick_param='off'):
    """
        Setting default panel tick parameters. Some of these are the defaults
        for matplotlib anyway, but specifying them for readability. Here are
        options and defaults for the parameters used below:

        axis : {'x', 'y', 'both'}; which axis to modify; default = 'both'
        which : {'major', 'minor', 'both'}; which ticks to modify;
                default = 'major'
        bottom, top, left, right : bool or {'on', 'off'}; ticks on or off;
        labelbottom, labeltop, labelleft, labelright : bool or {'on', 'off'}
     """

    # create the panel
    panel_rectangle = [panel_left, panel_bottom, panel_width, panel_height]
    panel = plt.axes(panel_rectangle)

    # Set tick parameters
    panel.tick_params(axis=axis_tick_param,
                      which=which_tick_param,
                      bottom=bottom_tick_param,
                      labelbottom=label_bottom_tick_param,
                      left=left_tick_param,
                      labelleft=label_left_tick_param,
                      right=right_tick_param,
                      labelright=label_right_tick_param,
                      top=top_tick_param,
                      labeltop=label_top_tick_param)

    return panel


def _generate_histogram_bin_patches(panel, bins, bin_values, horizontal=True):
    """This helper method generates the histogram that is added to the panel.

    In this case, horizontal = True applies to the mean quality histogram.
    So, horizontal = False only applies to the length histogram.
    """
    l_width = 0.0
    f_color = (0.5, 0.5, 0.5)
    e_color = (0, 0, 0)
    if horizontal:
        for step in np.arange(0, len(bin_values), 1):
            left = bins[step]
            bottom = 0
            width = bins[step + 1] - bins[step]
            height = bin_values[step]
            hist_rectangle = mplpatches.Rectangle((left, bottom), width, height,
                                                  linewidth=l_width,
                                                  facecolor=f_color,
                                                  edgecolor=e_color)
            panel.add_patch(hist_rectangle)
    else:
        for step in np.arange(0, len(bin_values), 1):
            left = 0
            bottom = bins[step]
            width = bin_values[step]
            height = bins[step + 1] - bins[step]

            hist_rectangle = mplpatches.Rectangle((left, bottom), width, height,
                                                  linewidth=l_width,
                                                  facecolor=f_color,
                                                  edgecolor=e_color)
            panel.add_patch(hist_rectangle)


def generate_histogram(panel, data_list, min_plot_val, max_plot_val,
                       bin_interval, hist_horizontal=True,
                       left_spine=True, bottom_spine=True,
                       top_spine=False, right_spine=False, x_label=None,
                       y_label=None):

    bins = np.arange(0, max_plot_val, bin_interval)

    bin_values, bins2 = np.histogram(data_list, bins)

    # hist_horizontal is used for quality
    if hist_horizontal:
        panel.set_xlim([min_plot_val, max_plot_val])
        panel.set_ylim([0, max(bin_values * 1.1)])
    # and hist_horizontal == Fale is for read length
    else:
        panel.set_xlim([0, max(bin_values * 1.1)])
        panel.set_ylim([min_plot_val, max_plot_val])

    # Generate histogram bin patches, depending on whether we're plotting
    # vertically or horizontally
    _generate_histogram_bin_patches(panel, bins, bin_values, hist_horizontal)

    panel.spines['left'].set_visible(left_spine)
    panel.spines['bottom'].set_visible(bottom_spine)
    panel.spines['top'].set_visible(top_spine)
    panel.spines['right'].set_visible(right_spine)

    if y_label is not None:
        panel.set_ylabel(y_label)
    if x_label is not None:
        panel.set_xlabel(x_label)

def generate_square_map(panel, data_frame, plot_min_y, plot_min_x,
                      plot_max_y, plot_max_x, color,
                      xcol, ycol, **kwargs):
    """This generates the heatmap panels using squares. Everything is
    quantized by ints.
    """
    panel.set_xlim([plot_min_x, plot_max_x])
    panel.set_ylim([plot_min_y, plot_max_y])
    tempdf = data_frame[[xcol, ycol]]
    data_frame = tempdf.astype(int)

    querystring = "{}<={} and {}<={}".format(plot_min_y, ycol, plot_min_x, xcol)
    print(" - Filtering squares with {}".format(querystring))
    square_this = data_frame.query(querystring)

    querystring = "{}<{} and {}<{}".format(ycol, plot_max_y, xcol, plot_max_x)
    print(" - Filtering squares with {}".format(querystring))
    square_this = square_this.query(querystring)

    counts = square_this.groupby([xcol, ycol]).size().reset_index(name='counts')
    for index, row in counts.iterrows():
        x_pos = row[xcol]
        y_pos = row[ycol]
        thiscolor = color(row["counts"]/(counts["counts"].max()))
        rectangle1=mplpatches.Rectangle((x_pos,y_pos),1,1,
                                        linewidth=0,\
                                        facecolor=thiscolor)
        panel.add_patch(rectangle1)

    all_counts = counts["counts"]
    return all_counts

def generate_heat_map(panel, data_frame, plot_min_y, plot_min_x,
                      plot_max_y, plot_max_x, color,
                      xcol, ycol, **kwargs):
    panel.set_xlim([plot_min_x, plot_max_x])
    panel.set_ylim([plot_min_y, plot_max_y])

    querystring = "{}<={} and {}<={}".format(plot_min_y, ycol, plot_min_x, xcol)
    print(" - Filtering hexmap with {}".format(querystring))
    hex_this = data_frame.query(querystring)

    querystring = "{}<{} and {}<{}".format(ycol, plot_max_y, xcol, plot_max_x)
    print(" - Filtering hexmap with {}".format(querystring))
    hex_this = hex_this.query(querystring)

    # This single line controls plotting the hex bins in the panel
    hex_vals = panel.hexbin(hex_this[xcol], hex_this[ycol], gridsize=49,
                            linewidths=0.0, cmap=color)
    for each in panel.spines:
        panel.spines[each].set_visible(False)

    counts = hex_vals.get_array()
    return counts

def generate_legend(panel, counts, color):
    # completely custom for more control
    panel.set_xlim([0, 1])
    panel.set_ylim([0, 1000])
    panel.set_yticks([int(x) for x in np.linspace(0, 1000, 6)])
    panel.set_yticklabels([int(x) for x in np.linspace(0, max(counts), 6)])
    for i in np.arange(0, 1001, 1):
        rgba = color(i / 1001)
        alpha = rgba[-1]
        facec = rgba[0:3]
        hist_rectangle = mplpatches.Rectangle((0, i), 1, 1,
                                              linewidth=0.0,
                                              facecolor=facec,
                                              edgecolor=(0, 0, 0),
                                              alpha=alpha)
        panel.add_patch(hist_rectangle)
    panel.spines['top'].set_visible(False)
    panel.spines['left'].set_visible(False)
    panel.spines['bottom'].set_visible(False)
    panel.yaxis.set_label_position("right")
    panel.set_ylabel('count')

def custommargin(df, **kwargs):
    rc.update_rcParams()

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

    # set the figure dimensions
    fig_width = 1.61 * 3
    fig_height = 1 * 3
    fig = plt.figure(figsize=(fig_width, fig_height))

    # set the panel dimensions
    heat_map_panel_width = fig_width * 0.5
    heat_map_panel_height = heat_map_panel_width * 0.62

    # find the margins to center the panel in figure
    fig_left_margin = fig_bottom_margin = (1 / 6)

    # lengthPanel
    y_panel_width = (1 / 8)

    # the color Bar parameters
    legend_panel_width = (1 / 24)

    # define padding
    h_padding = 0.02
    v_padding = 0.05

    # Set whether to include y-axes in histograms
    print(" - Setting panel options.", file = sys.stderr)
    if kwargs["Y_AXES"]:
        y_bottom_spine = True
        y_bottom_tick = 'on'
        y_bottom_label = 'on'
        x_left_spine = True
        x_left_tick = 'on'
        x_left_label = 'on'
        x_y_label = 'Count'
    else:
        y_bottom_spine = False
        y_bottom_tick = 'off'
        y_bottom_label = 'off'
        x_left_spine = False
        x_left_tick = 'off'
        x_left_label = 'off'
        x_y_label = None

    panels = []

    # Quality histogram panel
    print(" - Generating the x-axis panel.", file = sys.stderr)
    x_panel_left = fig_left_margin + y_panel_width + h_padding
    x_panel_width = heat_map_panel_width / fig_width
    x_panel_height = y_panel_width * fig_width / fig_height
    x_panel = generate_panel(x_panel_left,
                                fig_bottom_margin,
                                x_panel_width,
                                x_panel_height,
                                left_tick_param=x_left_tick,
                                label_left_tick_param=x_left_label)
    panels.append(x_panel)

    # y histogram panel
    print(" - Generating the y-axis panel.", file = sys.stderr)
    y_panel_bottom = fig_bottom_margin + x_panel_height + v_padding
    y_panel_height = heat_map_panel_height / fig_height
    y_panel = generate_panel(fig_left_margin,
                                  y_panel_bottom,
                                  y_panel_width,
                                  y_panel_height,
                                  bottom_tick_param=y_bottom_tick,
                                  label_bottom_tick_param=y_bottom_label)
    panels.append(y_panel)

    # Heat map panel
    heat_map_panel_left = fig_left_margin + y_panel_width + h_padding
    heat_map_panel_bottom = fig_bottom_margin + x_panel_height + v_padding
    print(" - Generating the heat map panel.", file = sys.stderr)
    heat_map_panel = generate_panel(heat_map_panel_left,
                                    heat_map_panel_bottom,
                                    heat_map_panel_width / fig_width,
                                    heat_map_panel_height / fig_height,
                                    bottom_tick_param='off',
                                    label_bottom_tick_param='off',
                                    left_tick_param='off',
                                    label_left_tick_param='off')
    panels.append(heat_map_panel)
    heat_map_panel.set_title(kwargs["title"])

    # Legend panel
    print(" - Generating the legend panel.", file = sys.stderr)
    legend_panel_left = fig_left_margin + y_panel_width + \
        heat_map_panel_width / fig_width + h_padding
    legend_panel_bottom = fig_bottom_margin + x_panel_height + v_padding
    legend_panel_height = heat_map_panel_height / fig_height
    legend_panel = generate_panel(legend_panel_left, legend_panel_bottom,
                                  legend_panel_width, legend_panel_height,
                                  bottom_tick_param='off',
                                  label_bottom_tick_param='off',
                                  left_tick_param='off',
                                  label_left_tick_param='off',
                                  right_tick_param='on',
                                  label_right_tick_param='on')
    panels.append(legend_panel)

    #
    # Everything above this is just to set up the panels
    #
    ##################################################################

    # Set max and min viewing window for the xaxis
    if kwargs["plot_max_x"]:
        plot_max_x = kwargs["plot_max_x"]
    else:
        if kwargs["square"]:
            plot_max_x = df[kwargs["xcol"]].max()
        plot_max_x = max(np.ceil(df[kwargs["xcol"]]))
    plot_min_x = kwargs["plot_min_x"]

    # Set x bin sizes
    if kwargs["xbin"]:
        x_bin_interval = kwargs["xbin"]
    else:
        # again, this is just based on what looks good from experience
        x_bin_interval = 1

    # Generate x histogram
    print(" - Generating the x-axis histogram.", file = sys.stderr)
    generate_histogram(panel = x_panel,
                       data_list = df[kwargs['xcol']],
                       min_plot_val = plot_min_x,
                       max_plot_val = plot_max_x,
                       bin_interval = x_bin_interval,
                       hist_horizontal = True,
                       x_label=kwargs['xcol'],
                       y_label=x_y_label,
                       left_spine=x_left_spine)

    # Set max and min viewing window for the y axis
    if kwargs["plot_max_y"]:
        plot_max_y = kwargs["plot_max_y"]
    else:
        if kwargs["square"]:
            plot_max_y = df[kwargs["ycol"]].max()
        else:
            plot_max_y = max(np.ceil(df[kwargs["ycol"]]))

    plot_min_y = kwargs["plot_min_y"]
    # Set y bin sizes
    if kwargs["ybin"]:
        y_bin_interval = kwargs["ybin"]
    else:
        y_bin_interval = 1

    # Generate y histogram
    print(" - Generating the y-axis histogram.", file = sys.stderr)
    generate_histogram(panel = y_panel,
                       data_list = df[kwargs['ycol']],
                       min_plot_val = plot_min_y,
                       max_plot_val = plot_max_y,
                       bin_interval = y_bin_interval,
                       hist_horizontal = False,
                       y_label = kwargs['ycol'],
                       bottom_spine = y_bottom_spine)

    # Generate heat map
    if kwargs["square"]:
        print(" - Generating the square heatmap.", file = sys.stderr)
        counts = generate_square_map(panel = heat_map_panel,
                                   data_frame = df,
                                   plot_min_y = plot_min_y,
                                   plot_min_x = plot_min_x,
                                   plot_max_y = plot_max_y,
                                   plot_max_x = plot_max_x,
                                   color = purple1,
                                   xcol = kwargs["xcol"],
                                   ycol = kwargs["ycol"])
    else:
        print(" - Generating the heatmap.", file = sys.stderr)
        counts = generate_heat_map(panel = heat_map_panel,
                                   data_frame = df,
                                   plot_min_y = plot_min_y,
                                   plot_min_x = plot_min_x,
                                   plot_max_y = plot_max_y,
                                   plot_max_x = plot_max_x,
                                   color = purple1,
                                   xcol = kwargs["xcol"],
                                   ycol = kwargs["ycol"])

    # Generate legend
    print(" - Generating the legend.", file = sys.stderr)
    generate_legend(legend_panel, counts, purple1)

    # inform the user of the plotting window if not quiet mode
    #if not kwargs["QUIET"]:
    #    print("""plotting in the following window:
    #    {0} <= Q-score (x-axis) <= {1}
    #    {2} <= length  (y-axis) <= {3}""".format(
    #        plot_min_x, plot_max_x, min_plot_val, max_plot_val),
    #        file=stderr)

    # Print image(s)
    if kwargs["output_base_name"] is None:
        file_base = "custommargin"
    else:
        file_base = kwargs["output_base_name"]

    print(" - Saving your images", file = sys.stderr)
    print_images(
        base =file_base,
        image_formats=kwargs["fileform"],
        dpi=kwargs["dpi"],
        no_timestamp = kwargs["no_timestamp"],
        transparent= kwargs["no_transparent"])

def run(args):
    print(args)
    if not opath.exists(args.input_file):
        raise IOError("The input file does not exist: {}".format(
            args.input_file))
    df = pd.read_csv(args.input_file, header='infer', sep='\t')
    # make sure that the column names that were specified are actually
    #  in the dataframe
    if args.xcol not in df.columns:
        raise IOError("""The x-column name that you specified, {}, is not in the
        dataframe column names: {}""".format(args.xcol, df.columns))
    if args.ycol not in df.columns:
        raise IOError("""The y-column name that you specified, {}, is not in the
        dataframe column names: {}""".format(args.ycol, df.columns))
    print(" - Successfully read csv file. Here are a few lines:",
          file = sys.stderr)
    print(df.head(), file = sys.stderr)
    print(" - Plotting {} on the x-axis".format(args.xcol),file=sys.stderr)
    print(df[args.xcol].head(), file = sys.stderr)
    print(" - Plotting {} on the y-axis".format(args.ycol),file=sys.stderr)
    print(df[args.ycol].head(), file = sys.stderr)
    custommargin(df=df.dropna(), **vars(args))
