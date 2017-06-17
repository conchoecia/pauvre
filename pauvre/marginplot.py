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
if platform.system() == 'Linux':
    matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import os.path as opath
from pauvre import rcparams
from pauvre.functions import parse_fastq_length_meanqual
from pauvre.stats import stats
from Bio import SeqIO

#logging
import logging
logger = logging.getLogger('pauvre')

def marginplot(args):
    length, meanQual = parse_fastq_length_meanqual(args.fastq)
    stats(args.fastq, length, meanQual)

    if args.maxlen:
        maxPlotLength = args.maxlen
    else:
        maxPlotLength = int(np.percentile(length, 99))

    if args.lengthbin:
        lengthBinInterval = args.lengthbin
    else:
        #Dividing by 80 is based on what looks good from experience
        lengthBinInterval = int(maxPlotLength/80)
    lengthBins=np.arange(0,maxPlotLength,lengthBinInterval)
    lengthBinsBig=np.arange(0,maxPlotLength,lengthBinInterval * 4)

    if args.maxqual:
        maxPlotQual = args.maxqual
    else:
        maxPlotQual = max(np.ceil(meanQual))
    if args.qualbin:
        qualBinInterval = args.qualbin
    else:
        #again, this is just based on what looks good from experience
        qualBinInterval = maxPlotQual/85
    qualBins=np.arange(0, maxPlotQual, qualBinInterval)
    qualBinsBig=np.arange(0, maxPlotQual, qualBinInterval*4)
    #250, 231, 34 light yellow
    #67, 1, 85
    #R=np.linspace(65/255,1,101)
    #G=np.linspace(0/255, 231/255, 101)
    #B=np.linspace(85/255, 34/255, 101)

    #R=65/255, G=0/255, B=85/255
    Rf=65/255
    Bf=85/255
    pdict = {'red':  ((0.0, Rf, Rf),
                       (1.0, Rf, Rf)),
              'green': ((0.0, 0.0, 0.0),
                        (1.0, 0.0, 0.0)),
              'blue':  ((0.0, Bf, Bf),
                        (1.0, Bf, Bf)),
              'alpha':  ((0.0, 0.0, 0.0),
                       (1.0, 1.0, 1.0))
             }
    # Now we will use this example to illustrate 3 ways of
    # handling custom colormaps.
    # First, the most direct and explicit:
    purple1 = LinearSegmentedColormap('Purple1', pdict)



    #make the pandas dataset to query
    df = pd.DataFrame(list(zip(length, meanQual)), columns=['length', 'meanQual'])
    #only keep the dataframes that are finite
    df = df.dropna()

    # I think I can delete this whole block
    # hist = {}
    # #create the histogram for plotting the interior
    # for lengthStep in np.arange(0,len(lengthBinsBig) - 1,1):
    #     for qualStep in np.arange(0, len(qualBinsBig) - 1 ,1):
    #         lengthMin = lengthBinsBig[lengthStep]
    #         lengthMax = lengthBinsBig[lengthStep + 1]
    #         qualMin = qualBinsBig[qualStep]
    #         qualMax = qualBinsBig[qualStep + 1]
    #         hist[(qualMin, lengthMin)] = len(df.query('length>={} and length<{} and meanQual>={} and meanQual<{}'.format(
    #             lengthMin, lengthMax, qualMin, qualMax)))

    #set the figure dimensions
    figWidth = 1.61*3
    figHeight = 1*3
    fig = plt.figure(figsize=(figWidth,figHeight))

    #set the panel dimensions
    panelWidth = figWidth * 0.5
    panelHeight = panelWidth/1.61

    #find the margins to center the panel in figure
    leftMargin = (1/6)*figWidth
    bottomMargin = (1/6)*figHeight

    #lengthPanel
    lengthPanelWidth = figWidth * (1/8)

    #the color Bar parameters
    colorWidth= figWidth * (1/24)

    x_lims = [0,10]
    y_lims = [0,2]

    panels = []

    qualPanel=plt.axes([(leftMargin + lengthPanelWidth)/figWidth, #left
                     bottomMargin/figHeight,    #bottom
                     panelWidth/figWidth,   #width
                     lengthPanelWidth/figHeight])     #height
    qualPanel.tick_params(axis='both',which='both',\
                       bottom='on', labelbottom='on',\
                       left='off', labelleft='off', \
                       right='off', labelright='off',\
                       top='off', labeltop='off')
    panels.append(qualPanel)

    lengthPanel=plt.axes([leftMargin/figWidth, #left
                     (bottomMargin + lengthPanelWidth)/figHeight,    #bottom
                     lengthPanelWidth/figWidth,   #width
                     panelHeight/figHeight])     #height
    lengthPanel.tick_params(axis='both',which='both',\
                       bottom='off', labelbottom='off',\
                       left='on', labelleft='on', \
                       right='off', labelright='off',\
                       top='off', labeltop='off')
    panels.append(lengthPanel)


    #set the left panel for 5' splice sites
    panel0=plt.axes([(leftMargin+lengthPanelWidth)/figWidth, #left
                     (bottomMargin + lengthPanelWidth)/figHeight,    #bottom
                     panelWidth/figWidth,   #width
                     panelHeight/figHeight])     #height
    panel0.tick_params(axis='both',which='both',\
                       bottom='off', labelbottom='off',\
                       left='off', labelleft='off', \
                       right='off', labelright='off',\
                       top='off', labeltop='off')
    panels.append(panel0)
    panel0.set_title(args.title)

    colorPanel=plt.axes([(leftMargin+ lengthPanelWidth + panelWidth + 0.1)/figWidth, #left
                     (bottomMargin + lengthPanelWidth)/figHeight,    #bottom
                     colorWidth/figWidth,   #width
                     (panelHeight )/figHeight])     #height
    colorPanel.tick_params(axis='both',which='both',\
                       bottom='off', labelbottom='off',\
                       left='off', labelleft='off', \
                       right='on', labelright='on',\
                       top='off', labeltop='off')
    panels.append(colorPanel)

    #plot the length histogram on y-axis
    lengthBinsValues,bins2=np.histogram(length,lengthBins)
    lengthPanel.set_ylim([0, maxPlotLength])
    lengthPanel.set_xlim([0, max(lengthBinsValues * 1.1)])

    for step in np.arange(0,len(lengthBinsValues),1):
        left=1
        bottom=lengthBins[step]
        width=lengthBinsValues[step]
        height=lengthBins[step+1]-lengthBins[step]

        rectangle1=mplpatches.Rectangle((left,bottom),width,height,\
                                        linewidth=0.0,\
                                        facecolor=(0.5,0.5,0.5),\
                                        edgecolor=(0,0,0))
        lengthPanel.add_patch(rectangle1)
    lengthPanel.spines['top'].set_visible(False)
    lengthPanel.spines['right'].set_visible(False)
    lengthPanel.spines['bottom'].set_visible(False)
    lengthPanel.set_ylabel('Read Length')

    #plot the length histogram on x-axis
    qualBinsValues,bins2=np.histogram(meanQual, qualBins)
    qualPanel.set_ylim([0, max(qualBinsValues * 1.1)])
    qualPanel.set_xlim([0, maxPlotQual])
    for step in np.arange(0,len(qualBinsValues),1):
        left=qualBins[step]
        bottom=0
        width=qualBins[step+1]-qualBins[step]
        height=qualBinsValues[step]

        rectangle1=mplpatches.Rectangle((left,bottom),width,height,\
                                        linewidth=0.0,\
                                        facecolor=(0.5,0.5,0.5),\
                                        edgecolor=(0,0,0))
        qualPanel.add_patch(rectangle1)
    qualPanel.spines['top'].set_visible(False)
    qualPanel.spines['right'].set_visible(False)
    qualPanel.spines['left'].set_visible(False)
    qualPanel.set_xlabel('Phred Quality')


    #plot the length histogram on x-axis
    hexThis = df.query('length<{} and meanQual<{}'.format(
                       maxPlotLength,maxPlotQual))
    #print(hexThis)

    panel0.set_xlim([0, maxPlotQual])
    panel0.set_ylim([0, maxPlotLength])
    #This single line controls plotting the hex bins in the panel
    hexvals = panel0.hexbin(hexThis['meanQual'], hexThis['length'], gridsize=49,
                            linewidths=0.0, cmap=purple1)
    counts = hexvals.get_array()

    # plot the colorbar
    # completely custom for more control
    colorPanel.set_xlim([0,1])
    colorPanel.set_ylim([0,1000])
    colorPanel.set_yticks([int(x) for x in np.linspace(0, 1000, 6)])
    colorPanel.set_yticklabels( [int(x) for x in np.linspace(0,max(counts), 6)] )
    for i in np.arange(0,1001,1):
        rgba = purple1(i/1001)
        alpha=rgba[-1]
        facec=rgba[0:3]
        rectangle1=mplpatches.Rectangle((0,i),1,1,\
                                    linewidth=0.0,\
                                    facecolor=facec,\
                                    edgecolor=(0,0,0),\
                                    alpha=alpha)
        colorPanel.add_patch(rectangle1)
    colorPanel.spines['top'].set_visible(False)
    colorPanel.spines['left'].set_visible(False)
    colorPanel.spines['bottom'].set_visible(False)
    colorPanel.yaxis.set_label_position("right")
    colorPanel.set_ylabel('Number of Reads')


    for each in panel0.spines:
        panel0.spines[each].set_visible(False)

    filebase = opath.splitext(opath.basename(args.fastq))[0]

    for each in args.fileform:
        outname = "{}.{}".format(filebase, each)
        if each == 'png':
            plt.savefig(outname, dpi = args.dpi, transparent=args.transparent)
        else:
            plt.savefig(outname, transparent=args.transparent)

def run(parser, args):
    marginplot(args)
