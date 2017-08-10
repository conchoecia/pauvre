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

"""This file contains things related to parsing and plotting GFF files"""

from matplotlib.path import Path
import matplotlib.patches as patches

global chevron_width
global arrow_width
arrow_width = 150
chevron_width = 75


def _plot_lff(left_df, right_df, colorMap, y_pos, bar_thickness):
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
    return patch

def _plot_rff(left_df, right_df, colorMap, y_pos, bar_thickness):
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
    return patch

def gffplot_horizontal(args, gff_object, track_width, start_y):
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
    print(plottable_features)
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
    while idone == False:
        print("im in the overlap-pairing while loop i={}".format(i))
        # look ahead at all of the elements that overlap with the ith element
        jdone = False
        j = 1
        this_set_minimum_index = i
        this_set_maximum_index = i
        while jdone == False:
            print("new i= {} j={} len={}".format(i, j, len_plottable))
            #print(len_plottable)
            #print(plottable_features)
            # first make sure that we haven't gone off the end of the dataframe
            # This is an edge case where i has a jth element that overlaps with it,
            #  and j is the last element in the plottable features.
            if i+j == len_plottable:
                # this checks for the case that i is the last element of the
                #  plottable features.
                # In both of the above cases, we are done with both the ith and
                #  the jth features.
                if i == len_plottable-1:
                    # this is the last analysis, so set idone to true
                    #  to finish after this
                    idone = True
                    # the last one can't be in its own group, so just add it solo
                    these_features = plottable_features.loc[this_set_minimum_index:this_set_maximum_index,].copy(deep=True)
                    plot_order.append(these_features.reset_index(drop=True))
                    break
                jdone = True
            else:
                # if the lmost of the next gene overlaps with the rmost of
                #  the current one, it overlaps and couple together
                if plottable_features.loc[i+j, 'lmost'] < plottable_features.loc[i, 'rmost']:
                    # note that this feature overlaps with the current
                    this_set_maximum_index = i+j
                    # ... and we need to look at the next in line
                    j += 1
                else:
                    i += 1 + (this_set_maximum_index - this_set_minimum_index)
                    #add all of the things that grouped together once we don't find any more groups
                    these_features = plottable_features.loc[this_set_minimum_index:this_set_maximum_index,].copy(deep=True)
                    plot_order.append(these_features.reset_index(drop=True))
                    jdone = True

    for feature_set in plot_order:
        # plot_feature_hori handles overlapping cases as well as normal cases
        patches = gffplot_feature_hori(feature_set, colorMap,
                                       start_y, bar_thickness)
        for each in patches:
               myPatches.append(each)
#    # Now we add all of the tRNAs to this to plot, do it last to overlay
#    #  everything else
#    tRNAs = gff_object.features.query("featType == 'tRNA'")
#    tRNAs.reset_index(inplace=True, drop = True)
#    print(tRNAs)
#    tRNA_bar_thickness = bar_thickness * (0.8)
#    tRNA_x_center  = x_center + ((tRNA_bar_thickness * (0.2))/2)
#    #tRNA_bar_thickness = bar_thickness
#    #tRNA_x_center  = x_center
#    print("sequence_length = {}".format(sequence_length))
#    angle_ranges = []
#    for i in range(0,len(tRNAs)):
#        this_feature = tRNAs.loc[i,]
#        min_angle = angleMap[min(this_feature.loc['start'],this_feature.loc['stop'])]
#        max_angle = angleMap[max(this_feature.loc['start'],this_feature.loc['stop'])]
#        angle_ranges.append((min_angle, max_angle))
#        patches = gffplot_feature_hori(this_feature, colorMap, tRNA_x_center,
#                               tRNA_bar_thickness, direction,
#                               pd.Series([]))
#        for each in patches:
#            myPatches.append(each)
#    print("angle ranges of tRNAs")
#    print(angle_ranges)
    return myPatches


def gffplot_feature_hori(feature_df, colorMap, y_pos, bar_thickness):
    """This plots the track for a feature, and if there is something for
    'this_feature_overlaps_feature', then there is special processing to
    add the white bar and the extra slope for the chevron
    """
    arrow_width = 80
    chevron_width = 40
    myPatches = []
    #if there is only one feature to plot, then just plot it
    if len(feature_df) == 1:
        print("plotting a single thing")
        #print(this_feature['name'], "is not overlapping")
        # This plots this shape:  ___________
        #                        |           \
        #                        |___________/
        if feature_df.loc[0,'strand'] == '+':
            direction = 1
        elif feature_df.loc[0, 'strand'] == '-':
            direction = -1

        verts = [(feature_df.loc[0, 'start'], y_pos + bar_thickness),
                 (feature_df.loc[0, 'stop'] - (arrow_width*direction), y_pos + bar_thickness),
                 (feature_df.loc[0, 'stop'], y_pos + (bar_thickness/2)),
                 (feature_df.loc[0, 'stop'] - (arrow_width*direction), y_pos),
                 (feature_df.loc[0, 'start'], y_pos),
                 (feature_df.loc[0, 'start'], y_pos + bar_thickness),
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
        for i in range(len(feature_df)):
            # this tests for which left type we're dealing with
            if i == 0:
                # type could be lff or lfr
                if feature_df.loc[i, 'strand'] == '+':
                    if feature_df.loc[i + 1, 'strand'] == '+':
                        # plot a lff type
                        myPatches.append(_plot_lff(feature_df.iloc[i,], feature_df.iloc[i+1,],
                                                   colorMap, y_pos, bar_thickness))
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
                        myPatches.append(_plot_rff(feature_df.iloc[i-1,], feature_df.iloc[i,],
                                                   colorMap, y_pos, bar_thickness))
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
    return myPatches
