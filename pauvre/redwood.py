#!/usr/bin/env python
# -*- coding: utf-8 -*-

# pauvre - just a pore PhD student's plotting package
# Copyright (c) 2016-2017 Darrin T. Schultz. All rights reserved.
# twitter @conchoecia
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

# SAM/BAM todo
# i/ for loop took 269 seconds ~ 4.48 minutes

#things to do
# - make each layer operate independently
#   - for each of these, make the program figure out the total length either from the GFF or from the bam file
# - make the error "Your query was too stringent and no reads resulted..." not give a traceback.
# - raise the proper error if the sam file has no header.
# - figure out how to update rcParams every time we run a program
# - figure out if the poretools logger is actually doing anything
# - Write a better docstring for how plotArc works.
# - Write a docstring for seqOrder method. I don't remember what it does
# - write a better function to get the alignment length of the sam/bam file
#    right now it opens the file twice and only gets the length the first time
# - drop columns by name, not by column number (samFile.drop...)
#   - Here's another that needs to be fixed: samFile.drop(samFile.columns[3], axis=1, inplace=True)
# - args
#   - get the filename from args
#   - set up a double-mapped mode to wrap reads around the 0th coordinate for
#      circular assemblies
#   - Make the r-dist something that the user can change.
# Getting Everything on the Same Track
#  - Make a list of features to plot `plot_order = []` or something like that
#  - First, go through the GFF features and come up with all of the things that
#     AREN'T tRNAs that overlap. Store each set of overlaps as its own set of
#     features. For things with no overlap, add those to the `plot_order` alone
#  - Iterate through all elements of `plot_order`, if all elements are forward
#     (start < stop), then draw the element at the end first with no modification,
#     then for every subsequent element draw a white border around the arrow.
#     The element same element should have the start of its arrow drawn to 1
#     degree before the start of the next one. This will make a chevron pattern
#     That will show both the start and stop of both reads accurately.
#     - If both elements are reverse, then do the same thing, but in reverse
#     - If one element is in reverse, but another is forward, then make all the
#       elements in the set half-width since there are too many possible
#       combinations to code reliably.

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rcParams
import platform
import itertools
import pandas as pd
import numpy as np
import pysam
import scipy as sp
import time
import os, sys

# import the pauvre rcParams
import pauvre.rcparams as rc
from pauvre.functions import GFFParse


# following this tutorial to install helvetica
# https://github.com/olgabot/sciencemeetproductivity.tumblr.com/blob/master/posts/2012/11/how-to-set-helvetica-as-the-default-sans-serif-font-in.md
global hfont
hfont = {'fontname':'Helvetica'}

#logging
import logging
logger = logging.getLogger('pauvre')

class BAMParse():
    """This class reads in a sam/bam file and constructs a pandas
    dataframe of all the relevant information for the reads to pass on
    and plot.
    """
    def __init__(self, filename, doubled):
        self.filename = filename
        self.doubled = doubled
        #determine if the file is bam or sam
        self.filetype = os.path.splitext(self.filename)[1]
        #throw an error if the file is not bam or sam
        if self.filetype not in ['.bam']:
            raise Exception("""You have provided a file with an extension other than
                            '.bam', please check your command-line arguments""")
        #now make sure there is an index file for the bam file
        if not os.path.exists("{}.bai".format(self.filename)):
            raise Exception("""Your .bam file is there, but it isn't indexed and
            there isn't a .bai file to go with it. Use
            'samtools index <yourfile>.bam' to fix it.""")
        #now open the file and just call it a sambam file
        filetype_dict = {'.sam': '', '.bam': 'b'}
        self.sambam = pysam.AlignmentFile(self.filename, "r{}".format(filetype_dict[self.filetype]))
        self.seqlength = self.sambam.lengths[0]
        self.true_seqlength = self.seqlength if not self.doubled else int(self.seqlength/2)
        self.raw_depthmap = self.get_depthmap()
        self.features = self.parse()
        self.features.sort_values(by=['POS','MAPLEN'], ascending=[True, False] ,inplace=True)
        self.features.reset_index(inplace=True)
        self.features.drop('index', 1, inplace=True)

    def get_depthmap(self):
        depthmap = [0] * self.true_seqlength
        for p in self.sambam.pileup():
            index = p.reference_pos
            if index >= self.true_seqlength:
                index -= self.true_seqlength
            depthmap[index] += p.nsegments
        return depthmap

    def parse(self):
        data = {'POS': [], 'MAPQ': [], 'TUPS': [] }
        for read in self.sambam.fetch():
           data['POS'].append(read.reference_start + 1)
           data['TUPS'].append(self.cigar_parse(read.cigartuples))
           data['MAPQ'].append(read.mapq)
        features = pd.DataFrame.from_dict(data, orient='columns')
        features['ALNLEN'] = features['TUPS'].apply(self.aln_len)
        features['TRULEN'] = features['TUPS'].apply(self.tru_len)
        features['MAPLEN'] = features['TUPS'].apply(self.map_len)
        features['POS'] =    features['POS'].apply(self.fix_pos)
        return features

    def cigar_parse(self, tuples):
        """
        arguments:
         <tuples> a CIGAR string tuple list in pysam format

        purpose:
         This function uses the pysam cigarstring tuples format and returns
         a list of tuples in the internal format, [(20, 'M'), (5, "I")], et
         cetera. The zeroth element of each tuple is the number of bases for the
         CIGAR string feature. The first element of each tuple is the CIGAR
         string feature type.

        There are several feature types in SAM/BAM files. See below:
         'M' - match
         'I' - insertion relative to reference
         'D' - deletion relative to reference
         'N' - skipped region from the reference
         'S' - soft clip, not aligned but still in sam file
         'H' - hard clip, not aligned and not in sam file
         'P' - padding (silent deletion from padded reference)
         '=' - sequence match
         'X' - sequence mismatch
         'B' - BAM_CBACK (I don't actually know what this is)

        """
        # I used the map values from http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment
        psam_to_char = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S',
                        5: 'H', 6: 'P', 7: '=', 8: 'X', 9: 'B'}
        return [(value, psam_to_char[feature]) for feature, value in tuples]

    def aln_len(self, TUPS):
        """
        arguments:
         <TUPS> a list of tuples output from the cigar_parse() function.

        purpose:
         This returns the alignment length of the read to the reference.
         Specifically, it sums the length of all of the matches and deletions.
         In effect, this number is length of the region of the reference sequence to
         which the read maps. This number is probably the most useful for selecting
         reads to visualize in the mapped read plot.
        """
        return sum([pair[0] for pair in TUPS if pair[1] not in ['S', 'H', 'I']])

    def map_len(self, TUPS):
        """
        arguments:
         <TUPS> a list of tuples output from the cigar_parse() function.

        purpose:
         This function returns the map length (all matches and deletions relative to
         the reference), plus the unmapped 5' and 3' hard/soft clipped sequences.
         This number is useful if you want to visualize how much 5' and 3' sequence
         of a read did not map to the reference. For example, poor quality 5' and 3'
         tails are common in Nanopore reads.
        """
        return sum([pair[0] for pair in TUPS if pair[1] not in ['I']])

    def tru_len(self, TUPS):
        """
        arguments:
         <TUPS> a list of tuples output from the cigar_parse() function.

        purpose:
         This function returns the total length of the read, including insertions,
         deletions, matches, soft clips, and hard clips. This is useful for
         comparing to the map length or alignment length to see what percentage of
         the read aligned to the reference.
        """
        return sum([pair[0] for pair in TUPS])

    def fix_pos(self, start_index):
        """
        arguments:
         <TUPS> a list of tuples output from the cigar_parse() function.

        purpose:
         When using a doubled SAMfile, any reads that start after the first copy
         of the reference risk running over the plotting window, causing the program
         to crash. This function corrects for this issue by changing the start site
         of the read.

        Note: this will probably break the program if not using a double alignment
        since no reads would map past half the length of the single reference
        """
        if self.doubled:
            if start_index > int(self.seqlength/2):
                return start_index - int(self.seqlength/2) - 1
            else:
                return start_index
        else:
            return start_index

def plotArc(start_angle, stop_angle, radius, width, **kwargs):
    """ write a docstring for this function"""
    numsegments = 100
    theta = np.radians(np.linspace(start_angle+90, stop_angle+90, numsegments))
    centerx = 0
    centery = 0
    x1 = -np.cos(theta) * (radius)
    y1 = np.sin(theta) * (radius)
    stack1 = np.column_stack([x1, y1])
    x2 = -np.cos(theta) * (radius + width)
    y2 = np.sin(theta) *  (radius + width)
    stack2 = np.column_stack([np.flip(x2, axis=0), np.flip(y2,axis=0)])
    #add the first values from the first set to close the polygon
    np.append(stack2, [[x1[0],y1[0]]], axis=0)
    arcArray = np.concatenate((stack1,stack2), axis=0)
    return patches.Polygon(arcArray, True, **kwargs), ((x1, y1), (x2, y2))

def fix_query_reflength(sequence_length, queries, doubled):
    """
    arguments:
     <sequence_length> This is the reference fasta length. It should be 2x the actual
               length of the reference since this program takes a sam file from
               a concatenated reference.
     <queries> This is a list of SQL-type query strings. This is generated
                from argparse.

    purpose:
     This function takes in a list of queries to use for read filtering
     for the redwood plot. It is often not advisable to plot all mapped reads
     since many of them are too small relative to the reference length. Also,
     the point of a death star plot is to show continuity of a circular
     reference, so short reads aren't very helpful there either.

     Currently, this function only recognizes the keyword argument 'reflength'.
    """
    if not doubled:
        sequence_length = int(sequence_length * 2)
    for i in range(len(queries)):
        if 'reflength' in queries[i].split():
            queries[i] = queries[i].replace('reflength', str(int(sequence_length/2)))

def cust_log(base, val):
    try:
        #val = np.log(val)/np.log(base)
        val = np.log(val)
    except:
        val = 0
    return val

def plot_histo(panel, args, angleMap, width_map, start_radius, track_radius,
               thisLog = False, double = False, ticks = False):
    """Plots a histogram given a width map. Width map must be the true length
    of the circular genome"""
    myPatches = []
    if double:
        #plot a line in the middle if doubled to distinguish the center
        mid_radius = start_radius+((track_radius - (track_radius/100))/2)
        arc, arcArray = plotArc(start_angle=0, stop_angle=360,
                          radius=mid_radius,
                          width=track_radius/100, fc='yellow')
        myPatches.append(arc)
        #this is only 1/100 the width of the band

    for i in range(len(width_map)):
        iStartIndex = i
        iStopIndex = i + 1
        iStartAngle = angleMap[iStartIndex]
        iStopAngle = angleMap[iStopIndex]
        if double:
            logwidth = -1 * track_radius * (np.log(width_map[i]-9)/np.log(max(width_map))) * 0.495
            log_start_radius = start_radius+((track_radius - (track_radius/100))/2)
            width = track_radius * (width_map[i]/max(width_map)) * 0.495
            normal_start_radius = start_radius+((track_radius - (track_radius/100))/2) + (track_radius/100)
            #first plot the log inside
            arc, arcArray = plotArc(start_angle=iStartAngle, stop_angle=iStopAngle,
                          radius=log_start_radius,
                          width=logwidth, fc='black')
            myPatches.append(arc)
            # now plot the normal outside
            arc, arcArray = plotArc(start_angle=iStartAngle, stop_angle=iStopAngle,
                          radius=normal_start_radius,
                          width=width, fc='black')
        else:
            if thisLog:
                base = 5
                width = track_radius * (cust_log(base, width_map[i]-2)/cust_log(base, max(width_map)))
                #alpha = cust_log(base, width_map[i])/cust_log(base, max(width_map))

            else:
                width = track_radius * (width_map[i]/max(width_map))
            arc, arcArray = plotArc(start_angle=iStartAngle, stop_angle=iStopAngle,
                          radius=start_radius,
                          width=width, fc='black')
        myPatches.append(arc)

    maxx = start_radius
    maxy = 0
    kerning = 12
    if ticks:
        # plot the scalebar
        maxval=max(width_map)
        xs = []
        ys = []
        xend = start_radius
        value_list = [10, 100, 1000]
        for value in value_list:
            centerAngle = 45
            this_radius = start_radius + (track_radius * (np.log(value-9)/np.log(maxval)))
            arc, arcArray = plotArc(start_angle=centerAngle-1, stop_angle=centerAngle+1,
                              radius=this_radius,
                              width=track_radius/25, fc='red')
            xs.append(arcArray[0][0][-1])
            ys.append(arcArray[0][1][-1])
            myPatches.append(arc)

        middle = np.mean(ys) - (kerning/2)
        new_ys = []
        if len(value_list) % 2 == 0:
            ystart = middle - (kerning/2) - ((len(new_ys) - 1) * kerning)
            yend = middle + (kerning/2) + ((len(new_ys) - 1) * kerning)
            for y in np.arange(ystart, yend+kerning, kerning):
                new_ys.append(y)
        else:
            ystart = middle - kerning
            yend = middle + kerning
            for y in np.arange(ystart, yend + kerning, kerning):
                new_ys.append(y)
        print("xs", xs)
        print("len xs", len(xs))
        print("values", value_list)
        print("new ys", new_ys)
        for i in range(len(value_list)):
            xlist = [xs[i], xend]
            ylist = [ys[i], new_ys[i]]
            panel.plot(xlist, ylist, ls='--', color='red', lw=0.75)
            panel.text(xend, new_ys[i],
                       value_list[i], fontsize = 12,
                       ha='left', va='center',
                       color = 'black', **hfont)
            panel.text(start_radius, new_ys[-1]+kerning,
                       "Read Depth ", fontsize = 10,
                       ha='center', va='center',
                       color = 'black', **hfont)


    return myPatches, start_radius + track_radius, panel

def plot_reads(args, angleMap, widthDict, samFiledf, start_radius,
               doubled = False, collapse = False, track_width=False,
               track_depth = False, thisLog = False):
    """this function plots the reads in a sam file.
    Outputs the patches, and the final_radius"""

    # there should be a different width dict with log scale
    #if args.log and thisLog:
    #    widthDict = {'M':0.8, # match
    #                 'I':0.95,  # insertion relative to reference
    #                 'D':0.05, # deletion relative to reference
    #                 'N':0.1,  # skipped region from the reference
    #                 'S':0.1,  # soft clip, not aligned but still in sam file
    #                 'H':0.1,  # hard clip, not aligned and not in sam file
    #                 'P':0.1,  # padding (silent deletion from padded reference)
    #                 '=':0.1,  # sequence match
    #                 'X':0.1}  # sequence mismatch

    append_radius = 0
    myPatches = []
    depth_map = []
    # if we need to collapse the reads so they aren't 1 read per line, keep track
    #  of plotted depth at a position with depth_map
    if collapse == True:
        if doubled:
            plotted_depth_map = [0] * int((len(angleMap)/2))
        else:
            plotted_depth_map = [0] * len(angleMap)

    # If we define a track width, then use the max read depth
    #  of the track to figure out the read_width
    read_width = 1
    if track_width:
        read_width = track_width/track_depth

    # This while loops cycles through the the list until everything is plotted
    #  We start at the first index of the file and count down. If we plot it, then
    #  we remove that read from the dataframe. Once we get to the last index,
    #  we attempt to plot and afterward reset the index of the dataframe and
    #  start from the beginning again.
    i = 0
    plotted = False
    #this is used to dtermine where we are during collapse
    current_row = 1
    if args.log and thisLog:
        current_row = 10
    skip = False
    # for printing out the progress of plotting,
    # get the original number of rows to make placeholders
    original_rownum_charcount = len(str(len(samFiledf)))
    direction = 'for'
    print(samFiledf)
    while len(samFiledf) > 0:
        stringTuples = samFiledf.loc[i, 'TUPS']
        mapLen = samFiledf.loc[i, 'MAPLEN']
        #print("i: {} and maplen: {}".format(i, mapLen))
        #Subtract one because BAM uses 1-based indexing but plotting uses 0.
        # I think I could avoid this in the future by changing the parse
        start_index = samFiledf.loc[i, 'POS'] - 1
        start_angle= angleMap[start_index]
        stop_index = 0
        stop_angle = angleMap[stop_index]
        if args.log and thisLog:
            log = 100
            read_width = track_width * ((cust_log(log, current_row + 1)/cust_log(100, track_depth + 1))\
                                        - (cust_log(log, current_row)/cust_log(100, track_depth + 1)))
            #print("log track depth: {}".format(np.log10(track_depth + 1)))
            #print("log current row: {}".format(np.log10(current_row + 1)))
        #now look in the plotted_depth_map to see if there is already a read
        # that overlaps with the current read we are considering.
        if collapse:
            for collapse_i in range(start_index, start_index+mapLen):
                #print("looking at {} and found {}".format(collapse_i, plotted_depth_map[collapse_i]))
                if plotted_depth_map[collapse_i] >= current_row:
                    skip = True
                    break
        #if we don't skip it, we're gonna plot it!
        if not skip:
            # Note that we are plotting something here
            if collapse:
                for collapse_i in range(start_index, start_index+mapLen):
                    plotted_depth_map[collapse_i] = plotted_depth_map[collapse_i] + 1
            for tup in stringTuples:
                if tup[1] == 'I':
                    #If there is an insertion, back up halfway and make plot the
                    # insertion to visually show a "bulge" with too much sequence.
                    # do not advance the start index to resume normal plotting
                    # after the insertion.
                    iStartIndex = start_index-int(tup[0]/2)
                    iStopIndex = iStartIndex + tup[0]
                    iStartAngle = angleMap[iStartIndex]
                    iStopAngle = angleMap[iStopIndex]
                    arc, arcArray =  plotArc(start_angle=iStartAngle, stop_angle=iStopAngle,
                                  radius=start_radius + append_radius,
                                  width=widthDict[tup[1]] * read_width, fc='black')
                else:
                    stop_index = start_index + tup[0]
                    stop_angle = angleMap[stop_index]
                    arc, arcArray =  plotArc(start_angle=start_angle, stop_angle=stop_angle,
                                  radius=start_radius + append_radius,
                                  width=widthDict[tup[1]] * read_width, fc='black')
                    start_index = stop_index
                    start_angle = angleMap[start_index]
                myPatches.append(arc)
            # If we're not collapsing the reads, we just advance one every row
            if not collapse:
                append_radius += read_width
            plotted = True
            #if we've ploted something, remove that read and reset the indices
            # since we will stay at the current index
            samFiledf.drop(i, inplace=True)
            samFiledf.reset_index(inplace=True)
            samFiledf.drop('index', 1, inplace=True)
            print("(countdown until done) rows: {0:0>{num}} dir: {1}\r".format(
                len(samFiledf), direction, num= original_rownum_charcount), end="")
        else:
            #if we weren't able to plot the current read, we will advance the
            # index and look for another read to plot
            # I tried to speed up the algorithm here by making the software 'smart'
            #  about looking forward for the next read, but it was ~50% slower
            i += 1
        skip = False
        # we will only reach this condition if we aren't collapsing, since all
        #  the reads will be removed before getting to this point
        if i >= len(samFiledf):
            i = 0
            #Once we've gone around a complete cycle, we can jump up to start
            # plotting the next row
            append_radius += read_width
            current_row += 1
            #every time we reset, reorganize so that we're now going in the
            # opposite direction to avoid skewing all the reads in the forward direction
            if args.interlace:
                if direction == 'for':
                    bav = {"by":['POS','MAPLEN'], "asc": [False, False]}
                    direction= 'rev'
                elif direction == 'rev':
                    bav = {"by":['POS','MAPLEN'], "asc": [True, False]}
                    direction = 'for'
                samFiledf.sort_values(by=bav["by"], ascending=bav['asc'],inplace=True)
                samFiledf.reset_index(inplace=True)
                samFiledf.drop('index', 1, inplace=True)
    print("\nfinal row is {}".format(current_row))
    return myPatches, start_radius + append_radius



def redwood(args):
    rc.update_rcParams()

    start = time.time()
    print(args)
    main_doubled = True if 'main' in args.doubled else False

    global sequence_length
    if args.main_bam:
        samFile = BAMParse(args.main_bam, main_doubled)
        sequence_length = samFile.seqlength 
        filename = samFile.filename
    else:
        sequence_length = GFFParse(args.gff).seqlen

    if sequence_length == 0:
        raise OSError("""You have used a SAM/BAM file with no header. Please add a header to
                 the file.""")
    # This stops numpy from printing numbers in scientific notation.
    np.set_printoptions(suppress=True)

    # this also needs to be changed depending on if it was a concatenated SAM
    # if doubled = true, then use linspace between 0,720
    # if doubled = false, then use linspace between 0, 360
    # on second thought, it might not be necessary to change this value even
    #  for doubled sequences
    global angleMap
    if 'main' in args.doubled:
        angleMap = np.linspace(0,720,sequence_length)
    else:
        angleMap = np.linspace(0,360,sequence_length)

    #these are the line width for the different cigar string flags.
    # usually, only M, I, D, S, and H appear in bwa mem output
    widthDict = {'M':0.45, # match
                 'I':0.9,  # insertion relative to reference
                 'D':0.05, # deletion relative to reference
                 'N':0.1,  # skipped region from the reference
                 'S':0.1,  # soft clip, not aligned but still in sam file
                 'H':0.1,  # hard clip, not aligned and not in sam file
                 'P':0.1,  # padding (silent deletion from padded reference)
                 '=':0.1,  # sequence match
                 'X':0.1}  # sequence mismatch
 
    ##################
    # make the redwood plot
    ##################

    figWidth = 5
    figHeight = 5

    # fig1 is an arbitrary name, only one figure currently
    fig_1 = plt.figure(figsize=(figWidth, figHeight))

    # This is the width and height of the plot in absolute
    #  values relative to the figWidth and figHeight.
    circleDiameter = 5.0

    # Center to plot
    leftMargin = (figWidth - circleDiameter)/2
    bottomMargin = leftMargin

    # There is one panel in this figure that contains
    #  the concentric circles that represent reads. In addition,
    #  the latest version of this program adds the annotation to the exterior
    panelCircle =  plt.axes([leftMargin/figWidth, #left
                             bottomMargin/figHeight,  #bottom
                             circleDiameter/figWidth, #width
                             circleDiameter/figHeight     #height
                             ],frameon=False)
    # Some of these are defaults and redundant, but are included
    #  for readability and convenience.
    panelCircle.tick_params(axis='both',which='both',\
                       bottom='off', labelbottom='off',\
                       left='off', labelleft='off', \
                       right='off', labelright='off',\
                       top='off', labeltop='off')

    # Simply the list in which patches will be stored.
    myPatches = []

    # Each read occupies a width of radius = 1. The max width occupied by any
    #  type of match to the read is 0.9 (Insertion). The line width of a match
    #  is 0.45. Look at the `widthDict` dictionary above to see the defined
    #  widths. I chose radius = 15 to start with because it looks decent in most
    #  scenarios, but maybe could use some tweaking in future iterations.
    start_radius_dict = {1: 10,
                         20: 15,
                         30: 20,
                         35: 90}
    if args.main_bam:
        print("The number of rows is: {}".format(len(samFile.features)))
        for key in sorted(start_radius_dict):
            if len(samFile.features) >= key:
                radius = start_radius_dict[key]
                print("Chose radius: {}".format(radius))
    else:
        radius = start_radius_dict[30]

    circle_fontsize = 10
    panelCircle.text(0, 0, "Position\n(bp)", fontsize = circle_fontsize,
                     ha='center', va='center',
                     color = 'black', **hfont)

    # now plot ticks around the interior, as well as the text
    if 'main' in args.doubled:
        real_length = int(sequence_length/2)
    else:
        real_length = sequence_length
    this_len = int(real_length / 1000) * 1000
    for i in range(0, this_len + 1000, 1000):
        startAngle = angleMap[i] - 0.75
        stopAngle = angleMap[i] + 0.75
        this_radius = radius * 0.93
        this_width =  radius * 0.04
        arc, arcArray = plotArc(start_angle=startAngle, stop_angle=stopAngle,
                      radius=this_radius,
                      width=this_width, fc='black')
        myPatches.append(arc)
        # now plot text if (val/1000) % 2 == 0
        if (i/1000) % 2 == 0:
            #the 0.98 gives some float
            x_pos = -np.cos(np.radians(angleMap[i]+90)) * this_radius
            y_pos =  np.sin(np.radians(angleMap[i]+90)) * this_radius
            rotation = 0
            if i == 0:
                y_pos =  np.sin(np.radians(angleMap[i]+90)) * this_radius * 0.98
                rotation = 0
                ha = 'center'
                va = 'top'
            if (angleMap[i] > 0 and angleMap[i] < 45):
                x_pos = -np.cos(np.radians(angleMap[i]+90 + 1.5)) * this_radius
                y_pos =  np.sin(np.radians(angleMap[i]+90 + 1.5)) * this_radius
                rotation = (90 - angleMap[i]) * 0.95
                ha = 'right'
                va = 'top'
            elif (angleMap[i] >= 45 and angleMap[i] < 67.5):
                rotation = 90 - angleMap[i]
                ha = 'right'
                va = 'top'
            elif (angleMap[i] >= 67.5 and angleMap[i] < 90):
                # subtracted an addtl degree because it looked bad otherwise
                x_pos = -np.cos(np.radians(angleMap[i]+90 - 1.75)) * this_radius
                y_pos =  np.sin(np.radians(angleMap[i]+90 - 1.75)) * this_radius
                rotation = 90 - angleMap[i]
                ha = 'right'
                va = 'top'
            elif (angleMap[i] >= 90 and angleMap[i] < 112.5):
                # added an addtl degree because it looked bad otherwise
                x_pos = -np.cos(np.radians(angleMap[i]+90 - 1.75)) * this_radius
                y_pos =  np.sin(np.radians(angleMap[i]+90 - 1.75)) * this_radius
                rotation = 90 - angleMap[i]
                ha = 'right'
                #bottom is not good
                va = 'center'
            elif (angleMap[i] >= 112.5 and angleMap[i] < 135):
                # added an addtl degree because it looked bad otherwise
                x_pos = -np.cos(np.radians(angleMap[i]+90 + 1.25)) * this_radius
                y_pos =  np.sin(np.radians(angleMap[i]+90 + 1.25)) * this_radius
                rotation = 90 - angleMap[i]
                ha = 'right'
                #bottom not good
                #center really not good
                va = 'bottom'
            elif (angleMap[i] >= 135 and angleMap[i] < 157.5):
                rotation = 90 - angleMap[i]
                ha = 'right'
                va = 'bottom'
            elif (angleMap[i] >= 157.5 and angleMap[i] < 180):
                x_pos = -np.cos(np.radians(angleMap[i]+90 - 1.0)) * this_radius
                y_pos =  np.sin(np.radians(angleMap[i]+90 - 1.0)) * this_radius
                rotation = 90 - angleMap[i]
                ha = 'right'
                va = 'bottom'
            elif (angleMap[i] >= 180 and angleMap[i] < 202.5):
                # added an addtl degree because it looked bad otherwise
                x_pos = -np.cos(np.radians(angleMap[i]+90 + 2.0)) * this_radius
                y_pos =  np.sin(np.radians(angleMap[i]+90 + 2.0)) * this_radius
                rotation = 270 - angleMap[i]
                ha = 'left'
                va = 'bottom'
            elif (angleMap[i] >= 202.5 and angleMap[i] < 225):
                x_pos = -np.cos(np.radians(angleMap[i]+90 + 2.25)) * this_radius
                y_pos =  np.sin(np.radians(angleMap[i]+90 + 2.25)) * this_radius
                rotation = 270 - angleMap[i]
                ha = 'left'
                va = 'bottom'
            elif (angleMap[i] >= 225 and angleMap[i] < 247.5):
                x_pos = -np.cos(np.radians(angleMap[i]+90 - 1.0)) * this_radius
                y_pos =  np.sin(np.radians(angleMap[i]+90 - 1.0)) * this_radius
                rotation = 270 - angleMap[i]
                ha = 'left'
                va = 'bottom'
            elif (angleMap[i] >= 247.5 and angleMap[i] < 270):
                x_pos = -np.cos(np.radians(angleMap[i]+90 - 3.0)) * this_radius
                y_pos =  np.sin(np.radians(angleMap[i]+90 - 3.0)) * this_radius
                rotation = 270 - angleMap[i]
                ha = 'left'
                va = 'bottom'
            elif (angleMap[i] >= 270 and angleMap[i] < 292.5):
                rotation = 270 - angleMap[i]
                ha = 'left'
                va = 'bottom'
            elif (angleMap[i] >= 292.5 and angleMap[i] < 315):
                x_pos = -np.cos(np.radians(angleMap[i]+90 + 0.25)) * this_radius
                y_pos =  np.sin(np.radians(angleMap[i]+90 + 0.25)) * this_radius
                rotation = 270 - angleMap[i]
                ha = 'left'
                va = 'top'
            elif (angleMap[i] >= 315 and angleMap[i] < 337.5):
                x_pos = -np.cos(np.radians(angleMap[i]+90 - 1.6)) * this_radius
                y_pos =  np.sin(np.radians(angleMap[i]+90 - 1.6)) * this_radius
                rotation = 270 - angleMap[i]
                ha = 'left'
                va = 'top'
            elif (angleMap[i] >= 337.5 and angleMap[i] < 360):
                x_pos = -np.cos(np.radians(angleMap[i]+90 - 5.0)) * this_radius
                y_pos =  np.sin(np.radians(angleMap[i]+90 - 5.0)) * this_radius
                rotation = 270 - angleMap[i]
                ha = 'left'
                va = 'top'
            print(" angleMap: {} value: {} rotation: {}".format(angleMap[i], i, rotation))
            text = "{}".format(i)
            if i == 0:
                text = "1/\n{}".format(real_length)
            panelCircle.text(x_pos, y_pos, text, fontsize = circle_fontsize,
                     ha=ha, va=va, color = 'black', rotation=rotation, **hfont)
            # Now add a legend


    # If the user wants to plot long reads, plot them
    if args.main_bam:
        #turn the query string into something usable,
        #  get rid of variables from argparse
        # Make sure we haven't told the program to not query
        if 'False' not in args.query:
            doubled = 'main' in args.doubled
            fix_query_reflength(sequence_length, args.query, doubled)
            # string all the queries together
            queryString = " and ".join(args.query)
            print("You are using this query string to filter reads:\n'{}'".format(queryString))
            samFile = samFile.features.query(queryString)
            if len(samFile) == 0:
                raise IOError("""Your query was too stringent and no reads resulted.  Please try
                again with a less stringent test.  Redwood plotter
                exiting""")

            # now determine how to sort the reads in the order they will be plotted
            if args.small_start == 'inside':
                ascend = True
            elif args.small_start == 'outside':
                ascend = False
                samFile.sort_values(by=args.sort, ascending=True, inplace=True)
                samFile = samFile.reset_index()
                samFile.drop('index', 1, inplace=True)

        #this plots the central rings from the sam file
        read_patches, radius = plot_reads(args, angleMap, widthDict, samFile,
                                          radius, doubled = True, collapse = False)
        myPatches = myPatches + read_patches

    # if the user would like to plot the annotation, plot it now. In the future,
    #  allow the user to select the order in which the individual elements are
    #  plotted. Since the annotation should have a fixed proportional size to
    #  the circle independent of the number of plotted reads, define a new
    #  radius context.

    if args.gff:
        print("the radius at the end of annotation is: {}".format(radius))
        panelCircle, gff_patches, gff_radius = plot_gff(args, panelCircle, args.gff, radius)
        myPatches = myPatches + gff_patches
        radius = gff_radius

    # it is helpful to be able to plot the RNAseq data along with the annotation.
    # plot that directly outside the annotation
    if args.rnaseq_bam:
        print("in RNAseq")
        rna_doubled = True if 'rnaseq' in args.doubled else False
        bamobject = BAMParse(args.rnaseq_bam, rna_doubled)
        samFile = bamobject.features
        track_width = radius * 0.15
        #read_patches, radius = plot_reads(args, angleMap,
        #                                  widthDict, samFile, radius,
        #                                  doubled = rna_doubled,
        #                                  collapse = True, track_width = track_width,
        #                                  track_depth = max(bamobject.raw_depthmap),
        #                                  thisLog = True)
        radius_orig=radius

        read_patches, radius, panelCircle = plot_histo(panelCircle, args, angleMap,
                                                       bamobject.get_depthmap(), radius,
                                                       track_width, thisLog = True, ticks = True)
        myPatches = myPatches + read_patches
        myPatches.append(arc)

    # The numseqs value is used to determine the viewing dimensions
    #  of the circle we will plot. It is scaled with the number of sequences
    #  that will be plotted. This value should probably be set last to
    #  accommodate other tracks, like annotation.
    min_radius = int(-5 - np.ceil(radius))
    max_radius = int(5 + np.ceil(radius))
    panelCircle.set_xlim([min_radius, max_radius])
    panelCircle.set_ylim([min_radius, max_radius])

    for patch in myPatches:
        panelCircle.add_patch(patch)

    end = time.time()
    print(end - start)
    #plt.show()
    for extension in args.fileform:
        plt.savefig('redwood.{}'.format(extension), dpi=args.dpi,
                    transparent=True, frameon=False)

def feature_set_direction(feature_set_df):
    """This function determines if the features in the dataframe passed here are
    all forward, all reverse, or mixed"""
    all_pos = all(feature_set_df['strand'] == '+')
    all_neg = all(feature_set_df['strand'] == '-')
    if all_pos:
        return '+'
    elif all_neg:
        return '-'
    else:
        return 'mixed'

def plot_feature(this_feature, colorMap, start_radius,
                 bar_thickness, direction, this_feature_overlaps_feature):
    """This plots the track for a feature, and if there is something for
    'this_feature_overlaps_feature', then there is special processing to
    add the white bar and the extra slope for the chevron
    """
    myPatches = []
    if this_feature_overlaps_feature.empty:
        iStartAngle = angleMap[this_feature['start']]
        iStopAngle = angleMap[this_feature['stop']] - 2
        arc, arcArray = plotArc(start_angle=iStartAngle,
                                  stop_angle=iStopAngle,
                                  radius = start_radius, width=bar_thickness,
                                  fc=colorMap[this_feature['featType']])
        myPatches.append(arc)
        #this bit plots the arrow triangles for the genes.
        #  Right now it makes each arrow only 1 degree in width and uses 100 segments to plot it.
        #  This resolution hasn't given me any artifacts thus far
        tStartAngle = angleMap[this_feature['stop']]-2
        tStopAngle = angleMap[this_feature['stop']]
        angles = np.linspace(tStartAngle,tStopAngle,100)
        widths = np.linspace(bar_thickness,0,100)
        for j in range(len(angles)-1):
            arc, arcArray =  plotArc(start_angle=tStartAngle,
                                      stop_angle=angles[j+1],
                                      radius=start_radius+(bar_thickness-widths[j+1])/2,
                                      width=widths[j+1],
                                      fc=colorMap[this_feature['featType']])
            myPatches.append(arc)
    else:
        # First, make a solid white bar
        iStartAngle = angleMap[this_feature['start']]
        iStopAngle = angleMap[this_feature_overlaps_feature['start']]
        arc, arcArray = plotArc(start_angle=iStartAngle,
                                  stop_angle=iStopAngle,
                                  radius = start_radius, width=bar_thickness,
                                  fc=colorMap['spacebar'])
        myPatches.append(arc)
        # Now, make the actual color bar
        iStartAngle = angleMap[this_feature['start']]
        iStopAngle = angleMap[this_feature_overlaps_feature['start']] - 1
        arc, arcArray =  plotArc(start_angle=iStartAngle,
                                  stop_angle=iStopAngle,
                                  radius = start_radius, width=bar_thickness,
                                  fc=colorMap[this_feature['featType']])
        myPatches.append(arc)
        #first plot a little pink bar for the outline
        tStartAngle = angleMap[this_feature_overlaps_feature['start']]
        tStopAngle = angleMap[this_feature['stop']]+1
        angles = np.linspace(tStartAngle,tStopAngle,100)
        widths = np.linspace(bar_thickness,0,100)
        for j in range(len(angles)-1):
            arc, arcArray = plotArc(start_angle=tStartAngle,
                                      stop_angle=angles[j+1],
                                      radius=start_radius+(bar_thickness-widths[j+1])/2,
                                      width=widths[j+1],
                                      fc=colorMap['spacebar'])
            myPatches.append(arc)
        #this bit plots the arrow triangles for the genes.
        #  Right now it makes each arrow only 1 degree in width and uses 100 segments to plot it.
        #  This resolution hasn't given me any artifacts thus far
        tStartAngle = angleMap[this_feature_overlaps_feature['start']]-1
        tStopAngle = angleMap[this_feature['stop']]
        angles = np.linspace(tStartAngle,tStopAngle,100)
        widths = np.linspace(bar_thickness,0,100)
        for j in range(len(angles)-1):
            arc, arcArray = plotArc(start_angle=tStartAngle,
                                    stop_angle=angles[j+1],
                                    radius=start_radius+(bar_thickness-widths[j+1])/2,
                                    width=widths[j+1],
                                    fc=colorMap[this_feature['featType']])
            myPatches.append(arc)
    return myPatches

def get_angles(name, center_angle, kerning_angle):
    num_chars = len(name.strip())
    if num_chars == 2:
        start_pos = center_angle
        stop_pos  = center_angle + kerning_angle
        return (start_pos, stop_pos)
    if num_chars % 2 == 0:
        start_pos = center_angle - (kerning_angle/2) - ((num_chars - 1)/2)
        stop_pos  = center_angle + (kerning_angle/2) + ((num_chars - 1)/2)
    else:
        start_pos = center_angle - (num_chars/2)
        stop_pos  = center_angle + (num_chars/2)
    angles = np.arange(start_pos, stop_pos + kerning_angle, kerning_angle)
    return angles

def plot_gff(args, panelCircle, gff_path, radius):
    #parse the gff file
    gffParser = GFFParse(gff_path)

    # Because this size should be relative to the circle that it is plotted next
    #  to, define the start_radius as the place to work from, and the width of
    #  each track.
    start_radius = radius
    track_width = radius * 0.15

    colorMap = {'gene': 'green', 'CDS': 'green', 'tRNA':'pink', 'rRNA':'red',
                'misc_feature':'purple', 'rep_origin':'orange', 'spacebar':'white',
                'ORF':'orange'}
    augment = 0
    bar_thickness = 0.9 * track_width
    # return these at the end
    myPatches=[]
    plot_order = []
    # this for loop relies on the gff features to already be sorted
    i = 0
    idone = False
    # we need to filter out the tRNAs since those are plotted last
    plottable_features = gffParser.features.query("featType != 'tRNA' and featType != 'region' and featType != 'source'")
    plottable_features.reset_index(inplace=True, drop=True)
    while idone == False:
        print("im in the overlap-pairing while loop i={}".format(i))
        # look ahead at all of the elements that overlap with the ith element
        jdone = False
        j = 1
        this_set_minimum_index = i
        this_set_maximum_index = i
        while jdone == False:
            print("new i= {} j={} len={}".format(i, j, len(plottable_features)))
            # first make sure that we haven't gone off the end of the dataframe
            if i+j == len(plottable_features):
                if i == len(plottable_features)-1:
                    # this is the last analysis, so set idone to true
                    #  to finish after this
                    idone = True
                    # the last one can't be in its own group, so just add it solo
                    these_features = plottable_features.loc[this_set_minimum_index:this_set_maximum_index,]
                    plot_order.append(these_features.reset_index(drop=True))
                    break
                jdone == True
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
                    these_features = plottable_features.loc[this_set_minimum_index:this_set_maximum_index,]
                    plot_order.append(these_features.reset_index(drop=True))
                    jdone = True

    for feature_set in plot_order:
        print(feature_set)
        direction = feature_set_direction(feature_set)
        print("direction = {}".format(direction))
        if direction == '+':
            for i in range(len(feature_set)-1, -1,-1):
                print("inside the plot for loop i = {}".format(i))
                this_feature = feature_set.loc[i,]
                print("got this single feature")
                #For the first element, just plot it normally
                if i == len(feature_set) - 1:
                    #plot the annotation
                    print("Im in the first plot thing")
                    patches = plot_feature(this_feature, colorMap, start_radius,
                                           bar_thickness, direction,
                                           pd.Series([]))
                    for each in patches:
                        myPatches.append(each)
                else:
                    print("now I'm plotting the other one")
                    overlapped_feature = feature_set.loc[i+1,]
                    print("this is the overlapped feature")
                    print(overlapped_feature)
                    patches = plot_feature(this_feature, colorMap, start_radius,
                                           bar_thickness, direction,
                                           overlapped_feature)
                    for each in patches:
                        myPatches.append(each)

            final_radius = start_radius + track_width
    # Now we add all of the tRNAs to this to plot, do it last to overlay
    #  everything else
    tRNAs = gffParser.features.query("featType == 'tRNA'")
    tRNAs.reset_index(inplace=True, drop = True)
    print(tRNAs)
    tRNA_bar_thickness = bar_thickness * (0.8)
    tRNA_start_radius  = start_radius + ((tRNA_bar_thickness * (0.2))/2)
    #tRNA_bar_thickness = bar_thickness
    #tRNA_start_radius  = start_radius
    print("sequence_length = {}".format(sequence_length))
    angle_ranges = []
    for i in range(0,len(tRNAs)):
        this_feature = tRNAs.loc[i,]
        min_angle = angleMap[min(this_feature.loc['start'],this_feature.loc['stop'])]
        max_angle = angleMap[max(this_feature.loc['start'],this_feature.loc['stop'])]
        angle_ranges.append((min_angle, max_angle))
        patches = plot_feature(this_feature, colorMap, tRNA_start_radius,
                               tRNA_bar_thickness, direction,
                               pd.Series([]))
        for each in patches:
            myPatches.append(each)
    print("angle ranges of tRNAs")
    print(angle_ranges)

    # now plot the text of the all the things. do this last
    #  to cover the previous things we plotted
    for i in range(len(gffParser.features.index)):
        if gffParser.features.loc[i, 'featType'] not in ['region', 'source', 'tRNA']:
            angles = []
            # to plot the centers, first figure out the center of the annotated
            #  gene by subtracting the end from the beginning
            name = gffParser.features.loc[i,'name']
            middle_position = int(gffParser.features.loc[i,'center'])
            # calculate the radius, because it changes depending on the track.
            #  Remember that the radius is the position at the bottom of the track,
            #  so we must use va=bottom.
            center_radius = radius + (gffParser.features.loc[i,'track'] * track_width) + (bar_thickness/2)
            # then, figure out the center angle of that position.
            #  I subtract one since I want to center this for the non-arrow
            #  part of the bar.
            center_angle = angleMap[middle_position]-2
            count = 0
            while count < 1:
                print()
                print(name)
                #figure out how many characters to plot
                kerning_angle = 2.5
                print("putting in {} into get angles with center angle = {}, kerning_angle = {}".format(name, center_angle, kerning_angle))
                angles = get_angles(name, center_angle, kerning_angle)
                char_min_angle = angles[0]-(kerning_angle / 2)
                char_max_angle = angles[-1] - (kerning_angle / 2)
                overlapping_angles = []
                # for every possible tRNA position, see if the minimum or
                #  max angle lies in the range of the text
                for tRNA_overlap_range in angle_ranges:
                    for tRNA_position in tRNA_overlap_range:
                        if (char_min_angle < tRNA_position and
                            tRNA_position < char_max_angle):
                            #overlapping_angles contains the tRNA overlapping angles now
                            overlapping_angles.append(tRNA_position)
                print("character angles {}".format(angles))
                print("overlapping angles {}".format(overlapping_angles))
                print("in while loop, center = {}".format(center_angle))
                if len(overlapping_angles) == 0:
                    count += 1
                else:
                    min_overlapping_angle = min(overlapping_angles)
                    print("min_overlapping_angle = {}".format(min_overlapping_angle))
                    max_overlapping_angle = max(overlapping_angles)
                    print("max_overlapping_angle = {}".format(max_overlapping_angle))
                    gene_start_angle = angleMap[gffParser.features.loc[i,'start']]
                    gene_stop_angle  = angleMap[gffParser.features.loc[i,'stop']]
                    start_dif = abs(min_overlapping_angle - gene_start_angle)
                    stop_dif =  abs(gene_stop_angle - max_overlapping_angle)
                    print("start_dif = {}, stop_dif = {}".format(start_dif, stop_dif))
                    if start_dif > stop_dif:
                        center_angle = min_overlapping_angle - ((min_overlapping_angle - gene_start_angle)/2)
                    else:
                        # the extra -2 is to compensate for the taper. This only
                        #  applies to the end of a gene stop. This code won't
                        #  work properly for - strand things currently
                        center_angle = gene_stop_angle - 2 - ((gene_stop_angle - max_overlapping_angle)/2)
                    count += 1
                    print("center angle after assignment = {}".format(center_angle))
                    #There is some overlap in the text and we need to
                    # figure out the min overlap and max overlap. Then figure
                    # out if the angle between min and start or max and end is greater
                    # use whichever is greater as the new place to put the text, then calculate the center
            angles = get_angles(name, center_angle, kerning_angle)
            text_angle = angles[-1] - angles[0]
            # this prints all the characters in their positions and rotates
            print_count = 0
            for j in range(len(angles)):
                this_width = gffParser.features.loc[i,'width']
                if not text_angle > (angleMap[this_width] - 1):
                    if center_angle > 95 and center_angle < 265:
                        rotation = (-1 * center_angle) - 180
                        this_char = name[len(angles) - 1 - j]
                    else:
                        rotation = -1 * center_angle
                        print("name: {} j: {}".format(name, j))
                        this_char = name[j]
                    this_angle = angles[j]
                    # now calculate the absolute x position
                    x_pos = -np.cos(np.radians(this_angle+90)) * center_radius
                    # now calculate the absolute y position
                    y_pos = np.sin(np.radians(this_angle+90)) * center_radius
                    #rotate the text if necessary
                    panelCircle.text(x_pos, y_pos, this_char, fontsize=10,
                                     ha='center', va='center',
                                     rotation = rotation,
                                     family = 'monospace',
                                     color = 'white')
                # this block handles the case where the text is too small to put
                #  parallel with the mitochondrial circle, so it is perpindicular
                else:
                    if print_count == 0:
                        if gffParser.features.loc[i,'strand'] == '+':
                            gene_start_angle = angleMap[gffParser.features.loc[i,'start']]
                            gene_stop_angle  = angleMap[gffParser.features.loc[i,'stop']]
                            center_angle = gene_stop_angle - ((gene_stop_angle-gene_start_angle)/2) 
                        # now calculate the absolute x position
                        x_pos = -np.cos(np.radians(center_angle+90)) * center_radius
                        # now calculate the absolute y position
                        y_pos = np.sin(np.radians(center_angle+90)) * center_radius
                        if (center_angle > 0 and center_angle < 180):
                            rotation = -1 * center_angle + 90
                        else:
                            rotation = -1 * center_angle - 90
                        panelCircle.text(x_pos, y_pos, name, fontsize=6,
                                     ha='center', va='center',
                                     rotation = rotation,
                                     family = 'monospace',
                                     color = 'white')

                        print_count += 1

    return panelCircle, myPatches, final_radius

def run(args):
    redwood(args)
