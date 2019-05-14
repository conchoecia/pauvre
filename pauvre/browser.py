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

# following this tutorial to install helvetica
# https://github.com/olgabot/sciencemeetproductivity.tumblr.com/blob/master/posts/2012/11/how-to-set-helvetica-as-the-default-sans-serif-font-in.md
global hfont
hfont = {'fontname':'Helvetica'}

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.patches as patches


import gffutils
import pandas as pd
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import numpy as np
import os
import pauvre.rcparams as rc
from pauvre.functions import GFFParse, print_images, timestamp
from pauvre import gfftools
from pauvre.lsi.lsi import intersection
from pauvre.bamparse import BAMParse
import progressbar
import platform
import sys
import time

# Biopython stuff
from Bio import SeqIO
import Bio.SubsMat.MatrixInfo as MI


class PlotCommand:
    def __init__(self, plotcmd, REF):
        self.ref = REF
        self.style_choices = []
        self.cmdtype = ""
        self.path = ""
        self.style = ""
        self.options = ""
        self._parse_cmd(plotcmd)

    def _parse_cmd(self, plotcmd):
        chunks = plotcmd.split(":")
        if chunks[0] == "ref":
            self.cmdtype = "ref"
            if len(chunks) < 2:
                self._len_error()
            self.path = self.ref
            self.style = chunks[1]
            self.style_choices = ["normal", "colorful"]
            self._check_style_choices()
            if len(chunks) > 2:
                self.options = chunks[2].split(",")
        elif chunks[0] in ["bam", "peptides"]:
            if len(chunks) < 3:
                self._len_error()
            self.cmdtype = chunks[0]
            self.path = os.path.abspath(os.path.expanduser(chunks[1]))
            self.style = chunks[2]
            if self.cmdtype == "bam":
                self.style_choices = ["depth", "reads"]
            else:
                self.style_choices = ["depth"]
            self._check_style_choices()
            if len(chunks) > 3:
                self.options = chunks[3].split(",")
        elif chunks[0] in ["gff3"]:
            if len(chunks) < 2:
                self._len_error()
            self.cmdtype = chunks[0]
            self.path = os.path.abspath(os.path.expanduser(chunks[1]))
            if len(chunks) > 2:
                self.options = chunks[2].split(",")


    def _len_error(self):
        raise IOError("""You selected {} to plot,
        but need to specify the style at least.""".format(self.cmdtype))
    def _check_style_choices(self):
        if self.style not in self.style_choices:
            raise IOError("""You selected {} style for
            ref. You must select from {}. """.format(
                self.style, self.style_choices))

global dna_color
dna_color = {"A": (81/255, 87/255, 251/255, 1),
             "T": (230/255, 228/255, 49/255, 1),
             "G": (28/255, 190/255, 32/255, 1),
             "C": (220/255, 10/255, 23/255, 1)}

#these are the line width for the different cigar string flags.
# usually, only M, I, D, S, and H appear in bwa mem output
global widthDict
widthDict = {'M':0.45, # match
             'I':0.9,  # insertion relative to reference
             'D':0.05, # deletion relative to reference
             'N':0.1,  # skipped region from the reference
             'S':0.1,  # soft clip, not aligned but still in sam file
             'H':0.1,  # hard clip, not aligned and not in sam file
             'P':0.1,  # padding (silent deletion from padded reference)
             '=':0.1,  # sequence match
             'X':0.1}  # sequence mismatch


global richgrey
richgrey = (60/255, 54/255, 69/255, 1)

def plot_ref(panel, chrid, start, stop, thiscmd):
    panel.set_xlim([start, stop])
    panel.set_ylim([-2.5, 2.5])
    panel.set_xticks([int(val) for val in np.linspace(start, stop, 6)])
    if thiscmd.style == "colorful":
        thisseq = ""
        for record in SeqIO.parse(thiscmd.ref, "fasta"):
            if record.id == chrid:
                thisseq = record.seq[start-1: stop]
        for i in range(len(thisseq)):
            left = start + i
            bottom = -0.5
            width = 1
            height = 1
            rect = patches.Rectangle((left, bottom),
                                     width, height,
                                     linewidth = 0,
                                     facecolor = dna_color[thisseq[i]] )
            panel.add_patch(rect)
    return panel

def safe_log10(value):
    try:
        logval = np.log10(value)
    except:
        logval = 0
    return logval

def plot_bam(panel, chrid, start, stop, thiscmd):
    bam = BAMParse(thiscmd.path)
    panel.set_xlim([start, stop])
    if thiscmd.style == "depth":
        maxdepth = max(bam.features_depthmap)
        maxdepthlog = safe_log10(maxdepth)
        if "log" in thiscmd.options:
            panel.set_ylim([-maxdepthlog, maxdepthlog])
            panel.set_yticks([int(val) for val in np.linspace(0, maxdepthlog, 2)])

        else:
            panel.set_yticks([int(val) for val in np.linspace(0, maxdepth, 2)])
            if "c" in thiscmd.options:
                panel.set_ylim([-maxdepth, maxdepth])
            else:
                panel.set_ylim([0, maxdepth])


        for i in range(len(bam.features_depthmap)):
            left = start + i
            width = 1
            if "c" in thiscmd.options and "log" in thiscmd.options:
                bottom = -1 * safe_log10(bam.features_depthmap[i])
                height = safe_log10(bam.features_depthmap[i]) * 2
            elif "c" in thiscmd.options and "log" not in thiscmd.options:
                bottom = -bam.features_depthmap[i]
                height = bam.features_depthmap[i] * 2
            else:
                bottom = 0
                height = bam.features_depthmap[i]
            if height > 0:
                rect = patches.Rectangle((left, bottom),
                                         width, height,
                                         linewidth = 0,
                                         facecolor = richgrey )
                panel.add_patch(rect)

    if thiscmd.style == "reads":
        #If we're plotting reads, we don't need y-axis
        panel.tick_params(bottom="off", labelbottom="off",
                          left ="off", labelleft = "off")
        reads = bam.features.copy()
        panel.set_xlim([start, stop])
        direction = "for"
        if direction == 'for':
            bav = {"by":['POS','MAPLEN'], "asc": [True, False]}
            direction= 'rev'
        elif direction == 'rev':
            bav = {"by":['POS','MAPLEN'], "asc": [True, False]}
            direction = 'for'
        reads.sort_values(by=bav["by"], ascending=bav['asc'],inplace=True)
        reads.reset_index(drop=True, inplace=True)

        depth_count = -1
        plotind = start
        while len(reads) > 0:
            #depth_count -= 1
            #print("len of reads is {}".format(len(reads)))
            potential = reads.query("POS >= {}".format(plotind))
            if len(potential) == 0:
                readsindex = 0
                #print("resetting plot ind from {} to {}".format(
                #    plotind, reads.loc[readsindex, "POS"]))
                depth_count -= 1

            else:
                readsindex = int(potential.index.values[0])
                #print("pos of potential is {}".format(reads.loc[readsindex, "POS"]))
            plotind = reads.loc[readsindex, "POS"]

            for TUP in reads.loc[readsindex, "TUPS"]:
                b_type = TUP[1]
                b_len = TUP[0]
                #plotting params
                # left same for all.
                left = plotind
                bottom = depth_count
                height = widthDict[b_type]
                width = b_len
                plot = True
                color = richgrey
                if b_type in ["H", "S"]:
                    """We don't plot hard or sort clips - like IGV"""
                    plot = False
                    pass
                elif b_type == "M":
                    """just plot matches normally"""
                    plotind += b_len
                elif b_type in ["D", "P", "=", "X"]:
                    """deletions get an especially thin line"""
                    plotind += b_len
                elif b_type == "I":
                    """insertions get a special purple bar"""
                    left = plotind - (b_len/2)
                    color = (200/255, 41/255, 226/255, 0.5)
                elif b_type == "N":
                    """skips for splice junctions, line in middle"""
                    bottom += (widthDict["M"]/2) - (widthDict["N"]/2)
                    plotind += b_len
                if plot:
                    rect = patches.Rectangle((left, bottom),
                                     width, height,
                                     linewidth = 0,
                                     facecolor = color )
                    panel.add_patch(rect)
            reads.drop([readsindex], inplace=True)
            reads.reset_index(drop = True, inplace=True)
        panel.set_ylim([depth_count, 0])

    return panel

def plot_gff3(panel, chrid, start, stop, thiscmd):

    db = gffutils.create_db(thiscmd.path, ":memory:")
    bottom = 0
    genes_to_plot = [thing.id
                     for thing in db.region(
                             region=(chrid, start, stop),
                             completely_within=False)
                     if thing.featuretype == "gene" ]
    #print("genes to plot are: " genes_to_plot)
    panel.set_xlim([start, stop])
    # we don't need labels on one of the axes
    #panel.tick_params(bottom="off", labelbottom="off",
    #                  left ="off", labelleft = "off")


    ticklabels = []
    for geneid in genes_to_plot:
        plotnow = False
        if "id" in thiscmd.options and geneid in thiscmd.options:
            plotnow = True
        elif "id" not in thiscmd.options:
            plotnow = True
        if plotnow:
            ticklabels.append(geneid)
            if db[geneid].strand == "+":
                panel = gfftools._plot_left_to_right_introns_top(panel, geneid, db,
                                                             bottom, text = None)
                bottom += 1
            else:
                raise IOError("""Plotting things on the reverse strand is
                not yet implemented""")
    #print("tick labels are", ticklabels)
    panel.set_ylim([0, len(ticklabels)])
    yticks_vals = [val for val in np.linspace(0.5, len(ticklabels) - 0.5, len(ticklabels))]
    panel.set_yticks(yticks_vals)
    print("bottom is: ", bottom)
    print("len tick labels is: ", len(ticklabels))
    print("intervals are: ", yticks_vals)
    panel.set_yticklabels(ticklabels)

    return panel

def browser(args):
    rc.update_rcParams()
    print(args)

    # if the user forgot to add a reference, they must add one
    if args.REF is None:
        raise IOError("You must specify the reference fasta file")

    # if the user forgot to add the start and stop,
    #  Print the id and the start/stop
    if args.CHR is None or args.START is None or args.STOP is None:
        print("""\n You have forgotten to specify the chromosome,
  the start coordinate, or the stop coordinate to plot.
  Try something like '-c chr1 --start 20 --stop 2000'.
  Here is a list of chromosome ids and their lengths
  from the provided reference. The minimum start coordinate
  is one and the maximum stop coordinate is the length of
  the chromosome.\n\nID\tLength""")
        for record in SeqIO.parse(args.REF, "fasta"):
            print("{}\t{}".format(record.id, len(record.seq)))
        sys.exit(0)

    if args.CMD is None:
        raise IOError("You must specify a plotting command.")

    # now we parse each set of commands
    commands = [PlotCommand(thiscmd, args.REF)
                for thiscmd in reversed(args.CMD)]

    # set the figure dimensions
    if args.ratio:
        figWidth = args.ratio[0] + 1
        figHeight = args.ratio[1] + 1
        #set the panel dimensions
        panelWidth = args.ratio[0]
        panelHeight =  args.ratio[1]

    else:
        figWidth = 7
        figHeight = len(commands) + 2
        #set the panel dimensions
        panelWidth = 5
        # panel margin x 2 + panel height = total vertical height
        panelHeight = 0.8
        panelMargin = 0.1

    figure = plt.figure(figsize=(figWidth,figHeight))

    #find the margins to center the panel in figure
    leftMargin = (figWidth - panelWidth)/2
    bottomMargin = ((figHeight - panelHeight)/2) + panelMargin

    plot_dict = {"ref": plot_ref,
                 "bam": plot_bam,
                 "gff3": plot_gff3
                 #"peptides": plot_peptides
                 }

    panels = []
    for i in range(len(commands)):
        thiscmd = commands[i]
        if thiscmd.cmdtype in ["gff3", "ref", "peptides"] \
           or thiscmd.style == "depth" \
           or "narrow" in thiscmd.options:
            temp_panelHeight = 0.5
        else:
            temp_panelHeight = panelHeight
        panels.append( plt.axes([leftMargin/figWidth, #left
                              bottomMargin/figHeight,    #bottom
                              panelWidth/figWidth,   #width
                              temp_panelHeight/figHeight])     #height
                      )
        panels[i].tick_params(axis='both',which='both',\
                           bottom='off', labelbottom='off',\
                           left='on', labelleft='on', \
                           right='off', labelright='off',\
                           top='off', labeltop='off')
        if thiscmd.cmdtype == "ref":
            panels[i].tick_params(bottom='on', labelbottom='on')



        #turn off some of the axes
        panels[i].spines["top"].set_visible(False)
        panels[i].spines["bottom"].set_visible(False)
        panels[i].spines["right"].set_visible(False)
        panels[i].spines["left"].set_visible(False)

        panels[i] = plot_dict[thiscmd.cmdtype](panels[i], args.CHR,
                                     args.START, args.STOP,
                                     thiscmd)

        bottomMargin = bottomMargin + temp_panelHeight + (2 * panelMargin)

    # Print image(s)
    if args.BASENAME is None:
        file_base = 'browser_{}.png'.format(timestamp())
    else:
        file_base = args.BASENAME
    path = None
    if args.path:
        path = args.path
    transparent = args.transparent
    print_images(
        base_output_name=file_base,
        image_formats=args.fileform,
        dpi=args.dpi,
        no_timestamp = kwargs["no_timestamp"],
        path = path,
        transparent=transparent)


def run(args):
    browser(args)
