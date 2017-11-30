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

# import the pauvre rcParams
import pandas as pd
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import numpy as np
import os
import pauvre.rcparams as rc
from pauvre.functions import GFFParse, print_images
from pauvre import gfftools
from pauvre.lsi.lsi import intersection
import progressbar
import platform
import time

# for the shuffling algorithm
from itertools import product

# Biopython stuff
from Bio import SeqIO
import Bio.SubsMat.MatrixInfo as MI

# following this tutorial to install helvetica
# https://github.com/olgabot/sciencemeetproductivity.tumblr.com/blob/master/posts/2012/11/how-to-set-helvetica-as-the-default-sans-serif-font-in.md
global hfont
hfont = {'fontname':'Helvetica'}

if platform.system() == 'Linux':
    matplotlib.use('agg')

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.patches as patches

def shuffle_optimize_gffs(args, GFFs):
    """This function takes in a list of GFF objects and reshuffles the
    individual files such that the resulting sequence of GFF files has
    the minimum number of intersections when plotting synteny

    if args.optimum_order, then the program will find the global minimum
    arrangement using the first GFF file as the anchor.

    if not args.optimum_order, then the program will find the local minimum
    shuffle between every input pair of GFF files to plot in the best way possible
    given the input order.

    returns a list of GFF files from which the user can calculate plotting coordinates
    """
    # we use the first-input gff as the topmost sequence,
    #  and then find the best synteny match for the remaining sequences
    shuffled_gffs = []
    if args.optimum_order:
        firstgff = GFFs[0]
        # save the first gff file unadultered
        shuffled_gffs.append(firstgff)
        nextgffs = GFFs[1:]
        while len(nextgffs) > 0:
            obs_list = []
            for i in range(len(nextgffs)):
                # every observation will be stored here as a tuple.
                #  zeroth element is the num intersections with the current gff
                #  first element is the index of nextgffs
                #  second element is the GFF object
                shuffles = nextgffs[i].shuffle()
                for k in range(len(shuffles)):
                    shuf = shuffles[k]
                    coords = firstgff.couple(shuf, this_y = 0, other_y = 1)
                    num_inters = len(intersection(coords))
                    obs_list.append((num_inters, i, shuf))
                    print(obs_list[-1])
            intersections, gffixs, shufs = zip(*obs_list)
            # get the index of the shuffled gff with the least number of
            #  intersections to the current one against which we are comparing
            print("intersections", intersections)
            selected_ix = intersections.index(min(intersections))
            # save this gff to shuffled gffs to use later for plotting
            shuffled_gffs.append(shufs[selected_ix])
            # remove the origin of the shuffled gff from nextgffs
            del nextgffs[gffixs[selected_ix]]
            print("global minimum was {} intersections".format(min(intersections)))
            # now update the firstgff to the latest shuffled one we collected
            firstgff = shufs[selected_ix]
    # plot the gff files in the order in which you input them,
    #  but shuffle them to the order with least intersections
    else:
        # first we need to find the best arrangement by finding the combinations
        #  that share the most unique genes
        genes_series = [GFFs[i].get_unique_genes() for i in range(len(GFFs))]
        combinations_indices = [0]
        remaining_indices = list(range(1, len(GFFs)))
        done = False
        biggest_intersection_index = -1
        biggest_intersection_value = 0
        current_remaining_indices_index = 0
        while not done:
            #get the len of the intersection
            #print("combinations_indices: {}".format(combinations_indices))
            #print("current_remaining_indices_index: {}".format(current_remaining_indices_index))
            #print("remaining_indices[current_remaining_indices_index]: {}".format(remaining_indices[current_remaining_indices_index]))
            #print("genes_series[remaining_indices[current_remaining_indices_index]]: {}".format(genes_series[remaining_indices[current_remaining_indices_index]]))
            this_intersection_value = len(genes_series[combinations_indices[-1]] &\
                                          genes_series[remaining_indices[current_remaining_indices_index]])
            if this_intersection_value > biggest_intersection_value:
                biggest_intersection_value = this_intersection_value
                biggest_intersection_index = current_remaining_indices_index
            if current_remaining_indices_index < len(remaining_indices)-1:
                current_remaining_indices_index += 1
            else:
                combinations_indices.append(remaining_indices[biggest_intersection_index])
                del remaining_indices[biggest_intersection_index]
                biggest_intersection_value = 0
                current_remaining_indices_index = 0
                biggest_intersection_index = -1
            if len(remaining_indices) == 0:
                done = True
        # The best order of genes with the most shared genes
        #I don't know if this is really that useful though since many species will overlap.
        # In a future implementation of this program it might be necessary to do sub-sorting of this list to get the lest number of line intersections
        print("The best gene combination is {}".format(combinations_indices))
        # now we rearrange the GFFs to the best order
        new_GFFs = [GFFs[i] for i in combinations_indices]
        # If we're adding another copy of the top one, add it here before shuffling
        if args.sandwich:
            new_GFFs.append(new_GFFs[0])
        shuffles = [new_GFFs[i].shuffle() for i in range(len(new_GFFs))]
        print([len(shuffles[i]) for i in range(len(shuffles))])
        cumulative_least_shuffled_value = 999999999999999999999999999999999999
        bar = progressbar.ProgressBar()
        for combination in bar(list(product(*shuffles))):
            num_intersections = []
            first_genes = [str(combination[i].features[combination[i].features['featType'].isin(['gene', 'rRNA', 'CDS', 'tRNA'])]['name'].head(n=1)).split()[1] for i in range(len(combination))]
            # skip to the next iteration if all the genes aren't the same
            if args.start_with_aligned_genes and len(set(first_genes)) != 1:
                continue
            for i in range(len(new_GFFs) - 1):
                j = i + 1
                #figure out the best shuffle the next sequence
                coords = combination[i].couple(combination[j], this_y = i, other_y = j)
                num_intersections.append(len(intersection(coords)))
            if sum(num_intersections) < cumulative_least_shuffled_value:
                shuffled_gffs = combination
                cumulative_least_shuffled_value = sum(num_intersections)
                print("\nnew least global intersections: {}".format(sum(num_intersections)))
    return shuffled_gffs

def black_colormap():
    zeroone = np.linspace(0, 1, 100)
    colorrange = [(0,0,0,x) for x in zeroone]
    minblosum = min(MI.blosum62.values())
    maxblosum = max(MI.blosum62.values())
    colormap = {i: colorrange[int(translate(i, minblosum, maxblosum, 0, 99))]
                for i in range(minblosum, maxblosum + 1, 1)}
    return colormap

def translate(value, left_min, left_max, right_min, right_max):
    """This code maps values from the left range and interpolates to the
    corresponding range on the right. This is used to translate the amino acid
    substition matrix scores to a scale between 0 and 1 for making alphamaps.

    I don't know if this works if the directionality of the ranges are swapped.
      IE [5, -10] mapped to [0, 1]

    args:
      <value> - the value in [<left_min>:<left_max>] to scale between
                 [<right_min>:<right_max>]
      <left_min> - the 'min' of the left (source) range
      <left_max> - the 'max' of the left (source) range
      <right_min> - the 'min' of the right (target) range
      <right_max>   the 'max' of the right (target) range

    output:
      the <value>(float) scaled between <right_min> and <right_max>
    """
    # Figure out how 'wide' each range is
    left_span = left_max - left_min
    right_span = right_max - right_min

    # Convert the left range into a 0-1 range (float)
    value_scaled = float(value - left_min) / float(left_span)

    # Convert the 0-1 range into a value in the right range.
    return right_min + (value_scaled * right_span)

def _samplename_warning(samplename, filename):
    raise Warning("""There is a sample in your fasta alignments that
         does not match the samplenames from the gff filenames. Please
         rename this samplename to not contain any spaces or underscores.
         IE for sample 'NC016', '>NC_016_-_ND6' will not work but
         '>NC016_-_ND6' will work.

         Erroneous name: {}
                   File: {}""".format(samplename, os.path.basename(filename)))

def _samplelength_warning(samplename, genename, featType, gfflen, alnlen):
    raise Warning("""The length of the protein alignment isn't the same as the
         length in the GFF file for the sample. Maybe you used a sequence in the
         alignment that is different from the annotation source? Check if the
         stop codons are deleted/inserted from either the GFF or alignment. The
         protein alignment length should be 3 less than the gff length if the
         stop codons were included in the gff annotation.

         Sample name: {}
           feat type: {}
           gene name: {}
          gff length: {}
          aln length: {}""".format(samplename, featType, genename, gfflen, alnlen))

def _nosample_warning(samplename, alngenename, gffnames):
    raise Warning("""One of the gff files doesn't contain a sequence that the
         alignment file indicates should be present. Either the alignment file
         is misnamed or the sequence name in the GFF file is not what you
         intended.

           Sample name: {}
         aln gene name: {}
             gff names: {}""".format(samplename, alngenename, gffnames))

def get_alignments(args):
    """
    this reads in all the alignments from the fasta directory.
    """
    # This is a dict object with key as
    filelist = {os.path.splitext(x)[0]:os.path.join(os.path.abspath(args.aln_dir), x)
                   for x in os.listdir(args.aln_dir)
                   if os.path.splitext(x)[1]}
    # one entry in seq_dict is:
    # {seqname: {"featType": featType,
    #            "seqs": {samplename: seq},
    #            "indices": {samplename: indices}}
    seqs_dict = {}
    # go through every gene in the genelist
    for genename in filelist:
        thisFeatType = ""
        seqs_list    = {}
        indices_list = {}
        for record in SeqIO.parse(filelist[genename], "fasta"):
            # get the sample name and make sure that the sample names match
            samplename = record.id.replace("_", " ").split()[0]
            if samplename not in args.samplenames:
                _samplename_warning(samplename, filelist[genename])
            # first, get the sample features
            samplegff = args.samplenames[samplename].features
            featType = samplegff.loc[samplegff['name'] == genename, 'featType'].to_string().split()[1]
            # now we determine if this is a prot alignment or a nucleotide aln
            if featType in ['gene', 'CDS']:
                final_seq = "".join([x*3 for x in record.seq])
            elif featType == 'rRNA':
                final_seq = str(record.seq)
            # we now need to verify that the protein sequence is
            #  the length of the gene in the gff file. Do this by removing
            #  gaps in the alignment
            gfffilt = samplegff.loc[samplegff['name'] == genename, 'width']
            if len(gfffilt) == 0:
                _nosample_warning(samplename, genename, list(samplegff['name']))
            gfflen = int(gfffilt)
            aln = final_seq.replace("-", "")
            alnlen = len(aln)
            if gfflen != alnlen:
                _samplelength_warning(samplename, genename, featType, gfflen, alnlen)
            # If we've made it this far without any errors, then incorporate the
            #  indices for each index
            #print("start_index", start_index)
            final_indices = [-1] * len(final_seq)
            # up until the next for loop, here we are determining which
            #  direction to move in. Reverse sequences decrease from the start
            strand = samplegff.loc[samplegff['name'] == genename, 'strand'].to_string().split()[1]
            if strand == '+':
                direction = 1
                start_index = int(samplegff.loc[samplegff['name'] == genename, 'start'])
            elif strand == '-':
                direction = -1
                start_index = int(samplegff.loc[samplegff['name'] == genename, 'stop'])
            for i in range(len(final_indices)):
                if final_seq[i] != '-':
                    final_indices[i] = start_index
                    start_index = start_index + (1 * direction)
            seqs_list[samplename] = final_seq
            if args.center_on:
                center_coord = int(args.samplenames[samplename].features.loc[args.samplenames[samplename].features['name'] == args.center_on, 'center'])
                indices_list[samplename] = np.array(final_indices) - center_coord
            else:
                indices_list[samplename] = final_indices
            thisFeatType = featType
        seqs_dict[genename] = {"featType": thisFeatType,
                              "seqs": seqs_list,
                               "indices": indices_list}
    return seqs_dict


def plot_synteny(seq1, ind1, seq2, ind2, y1, y2,
                 featType, matrix, cm, seqname):
    """This function plots all the lines for each"""
    myPatches = []
    colormap = {"COX1": '#c0d9ef',
                   "L": '#e8f1df',
                   "I": '#f7dedc',
                 "16S": '#ff2e00',
                 "12S": '#ffc239',
                 "cal": '#ffff54',
                "COX2": "#7fce66",
                 "ND2": "#00ae60",
                "COX3": "#00aeec",
                 "ND1": "#006fbb",
                   "*": "#ffffff",
                   "(": "#ded9c5",
                   "Q": "#ffc294",
                   "?": "#b5a2c4",
                 "ND4": "#968b5a",
                 "ND3": "#00fc65",
                "ND4L": "#00dcf0",
                 "ND6": "#ff994e",
                 "ND5": "#dc31e6",
                   "X": "#d8d8d8",
                   "G": "#abdce7",
                "CYTB": "#ff0059"}

    for i in range(len(seq1)):
        feat1 = seq1[i]
        feat2 = seq2[i]
        if feat1 != '-' and feat2 != '-':
            xs = []
            ys = []
            xs.append(ind1[i]) # top left
            ys.append(y1)
            xs.append(ind1[i] + 1) # top right
            ys.append(y1)
            xs.append(ind2[i] + 1) # bottom right
            ys.append(y2)
            xs.append(ind2[i]) #bottom left
            ys.append(y2)
            xs.append(ind1[i]) #top left
            ys.append(y1)
            alpha = 0.5
            if featType in ['CDS', 'gene']:
                try:
                    val = matrix[(feat1, feat2)]
                except:
                    val = matrix[(feat2, feat1)]
                color = cm[val]
                alpha = color[-1]
            elif featType == 'rRNA':
                if feat1 != feat2:
                    alpha=0
            color = colormap[seqname]
            stack1 = np.column_stack([xs, ys])
            myPatches.append(patches.Polygon(stack1, closed=True,
                                             color = color,
                                             alpha = alpha,
                                             lw=0))
    return myPatches

def synplot(args):
    rc.update_rcParams()
    print(args)
    GFFs = []
    for i in range(len(args.gff_paths)):
        gffpath = args.gff_paths[i]
        species = ""
        if args.gff_labels:
            species = args.gff_labels[i]
        GFFs.append(GFFParse(gffpath, args.stop_codons, species))

    # find the optimum shuffling pattern
    # and add a list of samplenames to the args
    optGFFs = shuffle_optimize_gffs(args, GFFs)
    # Make a sandwich for a circular comparison
    setattr(args, 'samplenames', {gff.samplename:gff for gff in optGFFs})

    # now get the cms and normalize
    #cms, normalize = gen_colormaps()
    cm = black_colormap()

    ## and we get the protein alignment objects
    # {seqname: {"featType": featType,
    #            "seqs": {samplename: seq},
    #            "indices": {samplename: indices}}
    seqs_dict = get_alignments(args)

    # now plot the lines as an example
    plt.style.use('BME163')

    # set the figure dimensions
    figWidth = 2.5*4
    figHeight = 5
    figure = plt.figure(figsize=(figWidth,figHeight))

    #set the panel dimensions
    panelWidth = 2.5 * 3
    panelHeight = 2.5

    #find the margins to center the panel in figure
    leftMargin = (figWidth - panelWidth)/2
    bottomMargin = ((figHeight - panelHeight)/2) + 0.25

    panel0=plt.axes([leftMargin/figWidth, #left
                     bottomMargin/figHeight,    #bottom
                     panelWidth/figWidth,   #width
                     panelHeight/figHeight])     #height
    panel0.tick_params(axis='both',which='both',\
                       bottom='on', labelbottom='off',\
                       left='off', labelleft='off', \
                       right='off', labelright='off',\
                       top='off', labeltop='off')
    #turn off some of the axes
    panel0.spines['top'].set_visible(False)
    panel0.spines['right'].set_visible(False)
    panel0.spines['left'].set_visible(False)

    # {seqname: {"featType": featType,
    #            "seqs": {samplename: seq},
    #            "indices": {samplename: indices}}
    allPatches = []
    for seqname in seqs_dict:
        #go through in order
        print(seqname)
        for i in range(0, len(optGFFs) - 1):
            samplei = optGFFs[i].samplename
            samplej = optGFFs[i+1].samplename
            if samplei in seqs_dict[seqname]["seqs"].keys() and\
               samplej in seqs_dict[seqname]["seqs"].keys():
                featType = seqs_dict[seqname]["featType"]
                seq1 = seqs_dict[seqname]["seqs"][samplei]
                ind1 = seqs_dict[seqname]["indices"][samplei]
                seq2 = seqs_dict[seqname]["seqs"][samplej]
                ind2 = seqs_dict[seqname]["indices"][samplej]
                # this is the top one, just leave it at the actual value since
                #  the base of the annotations start on the integer
                y1 = len(optGFFs) - 1 - i
                # this needs to be increased by the bar_thickness (0.9 * track_width in this case, or 0.09)
                y2 = len(optGFFs) - 2 - i
                myPatches = plot_synteny(seq1, ind1, seq2, ind2, y1, y2,
                                         featType, MI.blosum62, cm, seqname)
                for patch in myPatches:
                    allPatches.append(patch)

    print("len allPatches", len(allPatches))
    # this bit plots the simplified lines in the centers
    ## first we plot all the lines from the centers of matching genes.
    ##  This is temporary. Or maybe it should be a feature
    #verts = []
    #for i in range(len(optGFFs) - 1):
    #    j = i + 1
    #    coords = optGFFs[i].couple(optGFFs[j], this_y = len(optGFFs) - i, other_y = len(optGFFs) - i - 1)
    #    for coord in coords:
    #        verts.append(coord)

    #for vert in verts:
    #    xxyy = list(zip(*vert))
    #    panel0.plot(xxyy[0], xxyy[1])
    # now we plot horizontal lines showing the length of the mitochondrial sequence
    maxseqlen = 0
    for i in range(len(optGFFs)):
        gff = optGFFs[i]
        x_offset = 0
        if args.center_on:
            x_offset = -1 * int(gff.features.loc[gff.features['name'] == args.center_on, 'center'])
        panel0, patches = gfftools.gffplot_horizontal(figure, panel0, args, gff,
                                                      0.2, len(optGFFs) - i - 1 - (0.18/2),
                                                      x_offset = x_offset)
        seq_name = gff.features['sequence'].unique()[0]
        if args.gff_labels:
            seq_name = "$\it{{{0}}}$".format(gff.species)
        panel0.text(0 + x_offset, len(optGFFs) - i - 1 + (0.18/2),
                    seq_name, fontsize = 12,
                    ha='left', va='bottom',
                    color = 'black',
                    zorder = 100)

        if gff.seqlen > maxseqlen:
            maxseqlen = gff.seqlen
        xs = (1 + x_offset, gff.seqlen + x_offset)
        #ys = [len(optGFFs) - i - 1 + (0.09/2)]*2
        ys = [len(optGFFs) - i - 1]*2

        panel0.plot(xs, ys, color='black', zorder = -9)
        for patch in patches:
            allPatches.append(patch)

    for patch in allPatches:
        panel0.add_patch(patch)
    panel0.set_xlabel("position (bp)")

    #panel0.set_xlim([-15000, int(np.ceil(maxseqlen/1000)*1000)])
    #panel0.set_ylim([-0.5, 2.5])

    # This removes the text labels from the plot
    labels = [item.get_text() for item in panel0.get_xticklabels()]
    empty_string_labels = ['']*len(labels)
    panel0.set_xticklabels(empty_string_labels)

    # Print image(s)
    if args.BASENAME is None:
        file_base = 'synteny_{}.png'.format(timestamp())
    else:
        file_base = args.BASENAME
    transparent = args.TRANSPARENT
    print_images(file_base, args.fileform, args.dpi, transparent)


def timestamp():
    """
    Returns the current time in :samp:`YYYY-MM-DD HH:MM:SS` format.
    """
    return time.strftime("%Y%m%d_%H%M%S")


def run(args):
    synplot(args)
