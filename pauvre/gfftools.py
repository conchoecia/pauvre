#!/usr/bin/env python
# -*- coding: utf-8 -*-

# pauvre - a pore plotting package
# Copyright (c) 2016-2018 Darrin T. Schultz. All rights reserved.
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

"""This file contains things related to parsing and plotting GFF files"""

import copy
from matplotlib.path import Path
import matplotlib.patches as patches

global chevron_width
global arrow_width
global min_text
global text_cutoff

arrow_width = 80
chevron_width = 40
min_text = 550
text_cutoff = 150
import sys

global colorMap
colorMap = {'gene': 'green', 'CDS': 'green', 'tRNA':'pink', 'rRNA':'red',
                'misc_feature':'purple', 'rep_origin':'orange', 'spacebar':'white',
                'ORF':'orange'}

def _plot_left_to_right_introns(panel, geneid, db, y_pos, text = None):
    """ plots a left to right patch with introns when there are no intervening
    sequences to consider. Uses a gene id and gffutils database as input.
                         b
                 a    .-=^=-.    c
      1__________2---/   e   `---1__________2
      | #lff      \f            d| #lff      \
      | left to    \3            | left to    \3
      | right      /             | right      /
      5___________/4             5___________/4
    """
    #first we need to determine the number of exons
    bar_thickness = 0.75
    #now we can start plotting the exons
    exonlist = list(db.children(geneid, featuretype='CDS', order_by="start"))
    for i in range(len(exonlist)):
        cds_start = exonlist[i].start
        cds_stop =  exonlist[i].stop
        verts = [(cds_start, y_pos + bar_thickness), #1
                 (cds_stop - chevron_width, y_pos + bar_thickness), #2
                 (cds_stop, y_pos + (bar_thickness/2)), #3
                 (cds_stop - chevron_width, y_pos), #4
                 (cds_start, y_pos), #5
                 (cds_start, y_pos + bar_thickness), #1
        ]
        codes = [Path.MOVETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY,
        ]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, lw = 0,
                                  fc=colorMap['CDS'] )
        panel.add_patch(patch)

        # we must draw the splice junction
        if i < len(exonlist) - 1:
            next_start = exonlist[i+1].start
            next_stop =  exonlist[i+1].stop
            middle = cds_stop + ((next_start - cds_stop)/2)

            verts = [(cds_stop - chevron_width, y_pos + bar_thickness), #2/a
                     (middle, y_pos + 0.95), #b
                     (next_start, y_pos + bar_thickness), #c
                     (next_start, y_pos + bar_thickness - 0.05), #d
                     (middle, y_pos + 0.95 - 0.05), #e
                     (cds_stop - chevron_width, y_pos + bar_thickness -0.05), #f
                     (cds_stop - chevron_width, y_pos + bar_thickness), #2/a
                     ]
            codes = [Path.MOVETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.CLOSEPOLY,
                     ]
            path = Path(verts, codes)
            patch = patches.PathPatch(path, lw = 0,
                                      fc=colorMap['CDS'] )
            panel.add_patch(patch)

    return panel

def _plot_left_to_right_introns_top(panel, geneid, db, y_pos, text = None):
    """ slightly different from the above version such thatsplice junctions
    are more visually explicit.

    plots a left to right patch with introns when there are no intervening
    sequences to consider. Uses a gene id and gffutils database as input.
                            b
                    a    .-=^=-.    c
      1_____________2---/   e   `---1_____________2
      | #lff       /f              d| #lff       /
      | left to   /                 | left to   /
      | right    /                  | right    /
      4_________/3                  4_________/3
    """
    #first we need to determine the number of exons
    bar_thickness = 0.75
    #now we can start plotting the exons
    exonlist = list(db.children(geneid, featuretype='CDS', order_by="start"))
    for i in range(len(exonlist)):
        cds_start = exonlist[i].start
        cds_stop =  exonlist[i].stop
        verts = [(cds_start, y_pos + bar_thickness), #1
                 (cds_stop, y_pos + bar_thickness), #2
                 (cds_stop - chevron_width, y_pos), #4
                 (cds_start, y_pos), #5
                 (cds_start, y_pos + bar_thickness), #1
        ]
        codes = [Path.MOVETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY,
        ]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, lw = 0,
                                  fc=colorMap['CDS'] )
        panel.add_patch(patch)

        # we must draw the splice junction
        if i < len(exonlist) - 1:
            next_start = exonlist[i+1].start
            next_stop =  exonlist[i+1].stop
            middle = cds_stop + ((next_start - cds_stop)/2)

            verts = [(cds_stop-5, y_pos + bar_thickness), #2/a
                     (middle, y_pos + 0.95), #b
                     (next_start, y_pos + bar_thickness), #c
                     (next_start, y_pos + bar_thickness - 0.05), #d
                     (middle, y_pos + 0.95 - 0.05), #e
                     (cds_stop-5, y_pos + bar_thickness -0.05), #f
                     (cds_stop-5, y_pos + bar_thickness), #2/a
                     ]
            codes = [Path.MOVETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.CLOSEPOLY,
                     ]
            path = Path(verts, codes)
            patch = patches.PathPatch(path, lw = 0,
                                      fc=colorMap['CDS'] )
            panel.add_patch(patch)

    return panel

def _plot_lff(panel, left_df, right_df, colorMap, y_pos, bar_thickness, text):
    """ plots a lff patch
      1__________2      ____________
      | #lff      \     \ #rff      \
      | left for   \3     \ right for \
      | forward    /     / forward   /
      5___________/4    /___________/
    """
    #if there is only one feature to plot, then just plot it

    print("plotting lff")
    verts = [(left_df['start'], y_pos + bar_thickness), #1
             (right_df['start'] - chevron_width, y_pos + bar_thickness), #2
             (left_df['stop'], y_pos + (bar_thickness/2)), #3
             (right_df['start'] - chevron_width, y_pos), #4
             (left_df['start'], y_pos), #5
             (left_df['start'], y_pos + bar_thickness), #1
             ]
    codes = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY,
             ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, lw = 0,
                 fc=colorMap[left_df['featType']] )
    text_width = left_df['width']
    if text and text_width >= min_text:
        panel = _plot_label(panel, left_df, y_pos, bar_thickness)
    elif text and text_width < min_text and text_width >= text_cutoff:
        panel = _plot_label(panel, left_df,
                            y_pos, bar_thickness,
                            rotate = True, arrow = True)

    return panel, patch

def _plot_label(panel, df, y_pos, bar_thickness, rotate = False, arrow = False):
    # handles the case where a dataframe was passed
    fontsize = 8
    rotation = 0
    if rotate:
       fontsize = 5
       rotation = 90
    if len(df) == 1:
        x =((df.loc[0, 'stop'] - df.loc[0, 'start'])/2) + df.loc[0, 'start']
        y = y_pos + (bar_thickness/2)
        # if we need to center somewhere other than the arrow, need to adjust
        #  for the direction of the arrow
        # it doesn't look good if it shifts by the whole arrow width, so only
        #  shift by half the arrow width
        if arrow:
            if df.loc[0, 'strand'] == "+":
                shift_start = df.loc[0, 'start']
            else:
                shift_start = df.loc[0, 'start'] + (arrow_width/2)
            x =((df.loc[0, 'stop'] - (arrow_width/2) - df.loc[0, 'start'])/2) + shift_start
        panel.text(x, y,
                   df.loc[0, 'name'], fontsize = fontsize,
                   ha='center', va='center',
                   color = 'white', family = 'monospace',
                   zorder = 100, rotation = rotation)
    # and the case where a series was passed
    else:
        x = ((df['stop'] - df['start'])/2) + df['start']
        y = y_pos + (bar_thickness/2)
        if arrow:
            if df['strand'] == "+":
                shift_start = df['start']
            else:
                shift_start = df['start'] + (arrow_width/2)
            x =((df['stop'] - (arrow_width/2) - df['start'])/2) + shift_start
        panel.text(x, y,
                   df['name'], fontsize = fontsize,
                   ha='center', va='center',
                   color = 'white', family = 'monospace',
                   zorder = 100, rotation = rotation)

    return panel

def _plot_rff(panel, left_df, right_df, colorMap, y_pos, bar_thickness, text):
    """ plots a rff patch
      ____________      1__________2
      | #lff      \     \ #rff      \
      | left for   \    6\ right for \3
      | forward    /     / forward   /
      |___________/     /5__________/4
    """
    #if there is only one feature to plot, then just plot it

    print("plotting rff")
    verts = [(right_df['start'], y_pos + bar_thickness), #1
             (right_df['stop'] - arrow_width, y_pos + bar_thickness), #2
             (right_df['stop'], y_pos + (bar_thickness/2)), #3
             (right_df['stop'] - arrow_width, y_pos), #4
             (right_df['start'], y_pos), #5
             (left_df['stop'] + chevron_width, y_pos + (bar_thickness/2)), #6
             (right_df['start'], y_pos + bar_thickness), #1
             ]
    codes = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY,
             ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, lw = 0,
                 fc=colorMap[right_df['featType']] )
    text_width = right_df['width']
    if text and text_width >= min_text:
        panel = _plot_label(panel, right_df, y_pos, bar_thickness)
    elif text and text_width < min_text and text_width >= text_cutoff:
        panel = _plot_label(panel, right_df,
                            y_pos, bar_thickness, rotate = True)
    return panel, patch

def x_offset_gff(GFFParseobj, x_offset):
    """Takes in a gff object (a gff file parsed as a pandas dataframe),
    and an x_offset value and shifts the start, stop, center, lmost, and rmost.

    Returns a GFFParse object with the shifted values in GFFParse.features.
    """
    for columnname in ['start', 'stop', 'center', 'lmost', 'rmost']:
        GFFParseobj.features[columnname] = GFFParseobj.features[columnname] + x_offset
    return GFFParseobj

def gffplot_horizontal(figure, panel, args, gff_object,
                       track_width=0.2, start_y=0.1, **kwargs):
    """
    this plots horizontal things from gff files. it was probably written for synplot,
    as the browser does not use this at all.
    """
    # Because this size should be relative to the circle that it is plotted next
    #  to, define the start_radius as the place to work from, and the width of
    #  each track.
    colorMap = {'gene': 'green', 'CDS': 'green', 'tRNA':'pink', 'rRNA':'red',
                'misc_feature':'purple', 'rep_origin':'orange', 'spacebar':'white'}
    augment = 0
    bar_thickness = 0.9 * track_width
    # return these at the end
    myPatches=[]
    plot_order = []

    idone = False
    # we need to filter out the tRNAs since those are plotted last
    plottable_features = gff_object.features.query("featType != 'tRNA' and featType != 'region' and featType != 'source'")
    plottable_features.reset_index(inplace=True, drop=True)
    print(plottable_features)

    len_plottable = len(plottable_features)
    print('len plottable', len_plottable)
    # - this for loop relies on the gff features to already be sorted
    # - The algorithm for this loop works by starting at the 0th index of the
    #   plottable features (i).
    #   - It then looks to see if the next object (the jth) overlaps with the
    #     ith element.
    i = 0
    j = 1
    while i < len(plottable_features):
        if i + j == len(plottable_features):
            #we have run off of the df and need to include everything from i to the end
            these_features = plottable_features.loc[i::,].copy(deep=True)
            these_features = these_features.reset_index()
            print(these_features)
            plot_order.append(these_features)
            i = len(plottable_features)
            break
        print(" - i,j are currently: {},{}".format(i, j))
        stop = plottable_features.loc[i]["stop"]
        start = plottable_features.loc[i+j]["start"]
        print("stop: {}. start: {}.".format(stop, start))
        if plottable_features.loc[i]["stop"] <= plottable_features.loc[i+j]["start"]:
            print("    - putting elements {} through (including) {} together".format(i, i+j))
            these_features = plottable_features.loc[i:i+j-1,].copy(deep=True)
            these_features = these_features.reset_index()
            print(these_features)
            plot_order.append(these_features)
            i += 1
            j = 1
        else:
            j += 1

    #while idone == False:
    #    print("im in the overlap-pairing while loop i={}".format(i))
    #    # look ahead at all of the elements that overlap with the ith element
    #    jdone = False
    #    j = 1
    #    this_set_minimum_index = i
    #    this_set_maximum_index = i
    #    while jdone == False:
    #        print("new i= {} j={} len={}".format(i, j, len_plottable))
    #        print("len plottable in jdone: {}".format(len_plottable))
    #        print("plottable features in jdone:\n {}".format(plottable_features))
    #        # first make sure that we haven't gone off the end of the dataframe
    #        # This is an edge case where i has a jth element that overlaps with it,
    #        #  and j is the last element in the plottable features.
    #        if i+j == len_plottable:
    #            print("i+j == len_plottable")
    #            # this checks for the case that i is the last element of the
    #            #  plottable features.
    #            # In both of the above cases, we are done with both the ith and
    #            #  the jth features.
    #            if i == len_plottable-1:
    #                print("i == len_plottable-1")

    #                # this is the last analysis, so set idone to true
    #                #  to finish after this
    #                idone = True
    #                # the last one can't be in its own group, so just add it solo
    #                these_features = plottable_features.loc[this_set_minimum_index:this_set_maximum_index,].copy(deep=True)
    #                plot_order.append(these_features.reset_index(drop=True))
    #                break
    #            jdone = True
    #        else:
    #            print("i+j != len_plottable")
    #            # if the lmost of the next gene overlaps with the rmost of
    #            #  the current one, it overlaps and couple together
    #            if plottable_features.loc[i+j, 'lmost'] < plottable_features.loc[i, 'rmost']:
    #                print("lmost < rmost")
    #                # note that this feature overlaps with the current
    #                this_set_maximum_index = i+j
    #                # ... and we need to look at the next in line
    #                j += 1
    #            else:
    #                print("lmost !< rmost")
    #                i += 1 + (this_set_maximum_index - this_set_minimum_index)
    #                #add all of the things that grouped together once we don't find any more groups
    #                these_features = plottable_features.loc[this_set_minimum_index:this_set_maximum_index,].copy(deep=True)
    #                plot_order.append(these_features.reset_index(drop=True))
    #                jdone = True
    #        print("plot order is now: {}".format(plot_order))
    #        print("jdone: {}".format(str(jdone)))

    for feature_set in plot_order:
        # plot_feature_hori handles overlapping cases as well as normal cases
        panel, patches = gffplot_feature_hori(figure, panel, feature_set, colorMap,
                                              start_y, bar_thickness, text = True)
        for each in patches:
            print("there are {} patches after gffplot_feature_hori".format(len(patches)))
            print(each)
            myPatches.append(each)
            print("length of myPatches is: {}".format(len(myPatches)))

    # Now we add all of the tRNAs to this to plot, do it last to overlay
    #  everything else
    tRNAs = gff_object.features.query("featType == 'tRNA'")
    tRNAs.reset_index(inplace=True, drop = True)
    tRNA_bar_thickness = bar_thickness * (0.8)
    tRNA_start_y  = start_y + ((bar_thickness - tRNA_bar_thickness)/2)
    for i in range(0,len(tRNAs)):
        this_feature = tRNAs[i:i+1].copy(deep=True)
        this_feature.reset_index(inplace=True, drop = True)
        panel, patches = gffplot_feature_hori(figure, panel, this_feature, colorMap,
                             tRNA_start_y, tRNA_bar_thickness, text = True)
        for patch in patches:
            myPatches.append(patch)
    print("There are {} patches at the end of gffplot_horizontal()".format(len(myPatches)))
    return panel, myPatches

def gffplot_feature_hori(figure, panel, feature_df,
                         colorMap, y_pos, bar_thickness, text=True):
    """This plots the track for a feature, and if there is something for
    'this_feature_overlaps_feature', then there is special processing to
    add the white bar and the extra slope for the chevron
    """
    myPatches = []
    #if there is only one feature to plot, then just plot it
    if len(feature_df) == 1:
        #print("plotting a single thing: {} {}".format(str(feature_df['sequence']).split()[1],
        #                                              str(feature_df['featType']).split()[1] ))
        #print(this_feature['name'], "is not overlapping")
        # This plots this shape:  1_________2       2_________1
        #                        |  forward  \3   3/  reverse |
        #                        |5__________/4    \4________5|
        if feature_df.loc[0,'strand'] == '+':
            verts = [(feature_df.loc[0, 'start'], y_pos + bar_thickness),    #1
                     (feature_df.loc[0, 'stop'] - arrow_width, y_pos + bar_thickness), #2
                     (feature_df.loc[0, 'stop'], y_pos + (bar_thickness/2)), #3
                     (feature_df.loc[0, 'stop'] - arrow_width, y_pos),       #4
                     (feature_df.loc[0, 'start'], y_pos),                    #5
                     (feature_df.loc[0, 'start'], y_pos + bar_thickness)]    #1
        elif feature_df.loc[0,'strand'] == '-':
            verts = [(feature_df.loc[0, 'stop'], y_pos + bar_thickness),    #1
                     (feature_df.loc[0, 'start'] + arrow_width, y_pos + bar_thickness), #2
                     (feature_df.loc[0, 'start'], y_pos + (bar_thickness/2)), #3
                     (feature_df.loc[0, 'start'] + arrow_width, y_pos),       #4
                     (feature_df.loc[0, 'stop'], y_pos),                    #5
                     (feature_df.loc[0, 'stop'], y_pos + bar_thickness)]    #1
        feat_width = feature_df.loc[0,'width']
        if text and feat_width >= min_text:
            panel = _plot_label(panel, feature_df.loc[0,],
                                y_pos, bar_thickness)
        elif text and feat_width < min_text and feat_width >= text_cutoff:
            panel = _plot_label(panel, feature_df.loc[0,],
                                y_pos, bar_thickness,
                                rotate = True, arrow = True)

        codes = [Path.MOVETO,
                    Path.LINETO,
                    Path.LINETO,
                    Path.LINETO,
                    Path.LINETO,
                    Path.CLOSEPOLY]
        path = Path(verts, codes)
        print("normal path is: {}".format(path))
        # If the feature itself is smaller than the arrow, we need to take special measures to 
        if feature_df.loc[0,'width'] <= arrow_width:
            path = Path([verts[i] for i in [0,2,4,5]],
                        [codes[i] for i in [0,2,4,5]])
        patch = patches.PathPatch(path, lw = 0,
                  fc=colorMap[feature_df.loc[0, 'featType']] )
        myPatches.append(patch)
    # there are four possible scenarios if there are two overlapping sequences:
    #  ___________  ____________       ____________  ___________
    #  | #1       \ \ #1        \     / #2        / / #2        |
    #  | both seqs \ \ both seqs \   / both seqs / / both seqs  |
    #  | forward   / / forward   /   \ reverse   \ \  reverse   |
    #  |__________/ /___________/     \___________\ \___________|
    #  ___________  _____________      ____________ _ _________
    #  | #3       \ \ #3         |    / #2        _|  #2       \
    #  | one seq   \ \ one seq   |   / one seq   |_  one seq    \
    #  | forward    \ \ reverse  |   \ reverse    _|  forward   /
    #  |_____________\ \_________|    \__________|_ ___________/
    #
    #  These different scenarios can be thought of as different left/right
    #   flanking segment types.
    #   In the annotation #rff:
    #    - 'r' refers to the annotation type as being on the right
    #    - the first 'f' refers to the what element is to the left of this one.
    #      Since it is forward the 5' end of this annotation must be a chevron
    #    - the second 'f' refers to the right side of this element. Since it is
    #      forward it must be a normal arrow.
    #   being on the right
    #
    #   *LEFT TYPES*     *RIGHT TYPES*
    #  ____________      ____________
    #  | #lff      \     \ #rff      \
    #  | left for   \     \ right for \
    #  | forward    /     / forward   /
    #  |___________/     /___________/
    #  ___________      _____________
    #  | #lfr     \     \ #rfr       |
    #  | left for  \     \ right for |
    #  | reverse    \     \ reverse  |
    #  |_____________\     \_________|
    #    ____________    ___________
    #   / #lrr      /   / #rrr      |
    #  / left rev  /   / right rev  |
    #  \ reverse   \   \  reverse   |
    #   \___________\   \___________|
    #     ____________      __________
    #    / #lrf      _|   _| #rrf     \
    #   / left rev  |_   | _ right rev \
    #   \ forward    _|   _| forward   /
    #    \__________|    |____________/
    #
    # To properly plot these elements, we must go through each element of the
    #  feature_df to determine which patch type it is.
    elif len(feature_df) == 2:
        print("im in here feat len=2")
        for i in range(len(feature_df)):
            # this tests for which left type we're dealing with
            if i == 0:
                # type could be lff or lfr
                if feature_df.loc[i, 'strand'] == '+':
                    if feature_df.loc[i + 1, 'strand'] == '+':
                        # plot a lff type
                        panel, patch = _plot_lff(panel, feature_df.iloc[i,], feature_df.iloc[i+1,],
                                                 colorMap, y_pos, bar_thickness, text)
                        myPatches.append(patch)
                    elif feature_df.loc[i + 1, 'strand'] == '-':
                        #plot a lfr type
                        raise IOError("can't plot {} patches yet".format("lfr"))
                 # or type could be lrr or lrf
                elif feature_df.loc[i, 'strand'] == '-':
                    if feature_df.loc[i + 1, 'strand'] == '+':
                        # plot a lrf type
                        raise IOError("can't plot {} patches yet".format("lrf"))
                    elif feature_df.loc[i + 1, 'strand'] == '-':
                        #plot a lrr type
                        raise IOError("can't plot {} patches yet".format("lrr"))
            # in this case we're only dealing with 'right type' patches
            elif i == len(feature_df) - 1:
                # type could be rff or rfr
                if feature_df.loc[i-1, 'strand'] == '+':
                    if feature_df.loc[i, 'strand'] == '+':
                        # plot a rff type
                        panel, patch = _plot_rff(panel, feature_df.iloc[i-1,], feature_df.iloc[i,],
                                                 colorMap, y_pos, bar_thickness, text)
                        myPatches.append(patch)
                    elif feature_df.loc[i, 'strand'] == '-':
                        #plot a rfr type
                        raise IOError("can't plot {} patches yet".format("rfr"))
                # or type could be rrr or rrf
                elif feature_df.loc[i-1, 'strand'] == '-':
                    if feature_df.loc[i, 'strand'] == '+':
                        # plot a rrf type
                        raise IOError("can't plot {} patches yet".format("rrf"))
                    elif feature_df.loc[i, 'strand'] == '-':
                        #plot a rrr type
                        raise IOError("can't plot {} patches yet".format("rrr"))
    return panel, myPatches
