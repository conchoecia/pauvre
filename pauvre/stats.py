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

# TODO: make a function to nicely print out the pandas dataframes

"""
Program: pauvre stats

example usage/output:
 - For example, if I want to know the number of bases and reads with AT LEAST
   PHRED score 5 and AT LEAST length of 500, I can run the program as you
   see below and look at the cells highlighted with <braces>.
 - We find that there are 968,653 basepairs contained in 859 reads that
   fit meanReadQuality >= Q5, readLen >= 500.
 > pauvre marginplot --fastq miniDSMN15.fastq

 >numReads: 1000
 >numBasepairs: 1029114
 >meanLen: 1029.114
 >medianLen: 875.5
 >minLen: 11
 >maxLen: 5337
 >N50: 1278
 >L50: 296

 >                      Basepairs >= bin by mean PHRED and length
 >minLen       Q0       Q5     Q10     Q15   Q17.5    Q20  Q21.5   Q25  Q25.5  Q30
 >     0  1029114  1010681  935366  429279  143948  25139   3668  2938   2000    0
 >   500   984212  <968653> 904787  421307  142003  24417   3668  2938   2000    0
 >  1000   659842   649319  616788  300948  103122  17251   2000  2000   2000    0
 > et cetera...
 >              Number of reads >= bin by mean Phred+Len
 >minLen    Q0   Q5  Q10  Q15  Q17.5  Q20  Q21.5  Q25  Q25.5  Q30
 >     0  1000  969  865  366    118   22      3    2      1    0
 >   500   873 <859> 789  347    113   20      3    2      1    0
 >  1000   424  418  396  187     62   11      1    1      1    0
 > et cetera...
"""


from pauvre.functions import parse_fastq_length_meanqual
import os
import pandas as pd
import numpy as np


def stats(df, fastqName, histogram=False):
    """
    arguments:
     <df>
      - a pandas dataframe with cols ['length', 'meanQual'] that contains read
         length and quality information
     <fastqName>
      - just the name of the fastq file. This is used for printing out headers
         for the data tables.

    purpose:
     Calculate and output various stats about this fastq file. Currently
     this program reports:
       - Total number of reads
       - Total number of basepairs
       - mean
       - median
       - min
       - max
       - N50
       - A table of num basepairs filtered by min mean PHRED and length
       - A table of num reads filtered by the same as above ^

    example usage/output:
     - For example, if I want to know the number of bases and reads with AT LEAST
       PHRED score 5 and AT LEAST length of 500, I can run the program as you
       see below and look at the cells highlighted with <braces>.
     - We find that there are 968,653 basepairs contained in 859 reads that
       fit meanReadQuality >= Q5, readLen >= 500.
     > pauvre marginplot --fastq miniDSMN15.fastq

     >numReads: 1000
     >numBasepairs: 1029114
     >meanLen: 1029.114
     >medianLen: 875.5
     >minLen: 11
     >maxLen: 5337
     >N50: 1278
     >L50: 296

     >                      Basepairs >= bin by mean PHRED and length
     >minLen       Q0       Q5     Q10     Q15   Q17.5    Q20  Q21.5   Q25  Q25.5  Q30
     >     0  1029114  1010681  935366  429279  143948  25139   3668  2938   2000    0
     >   500   984212  <968653> 904787  421307  142003  24417   3668  2938   2000    0
     >  1000   659842   649319  616788  300948  103122  17251   2000  2000   2000    0
     > et cetera...
     >              Number of reads >= bin by mean Phred+Len
     >minLen    Q0   Q5  Q10  Q15  Q17.5  Q20  Q21.5  Q25  Q25.5  Q30
     >     0  1000  969  865  366    118   22      3    2      1    0
     >   500   873 <859> 789  347    113   20      3    2      1    0
     >  1000   424  418  396  187     62   11      1    1      1    0
     > et cetera...
    """
    fastqBase = os.path.basename(fastqName)

    analysis_sizes = [0, 1000, 5000, 10000]
    print_string = ""
    for this_size in analysis_sizes:
        these_lengths = df.loc[df["length"] >= this_size, "length"]
        if len(these_lengths) > 0:
            print_string += "# Fastq stats for {}, reads >= {}bp\n".format(fastqBase, this_size)
            print_string += "numReads: {}\n".format(len(these_lengths))
            print_string += "%totalNumReads: {0:.2f}\n".format((len(these_lengths) / len(df)) * 100)
            print_string += "numBasepairs: {}\n".format(sum(these_lengths))
            print_string += "%totalBasepairs: {0:.2f}\n".format(
                (sum(these_lengths) / sum(df["length"])) * 100)
            print_string += "meanLen: {}\n".format(np.mean(these_lengths))
            print_string += "medianLen: {}\n".format(np.median(these_lengths))
            print_string += "minLen: {}\n".format(min(these_lengths))
            maxLen = max(these_lengths)
            print_string += "maxLen: {}\n".format(maxLen)

            # calculate the N50
            fiftypercent = 0.5 * sum(these_lengths)
            lenSum = 0
            count = 0
            for val in sorted(these_lengths, reverse=True):
                lenSum += val
                count += 1
                if lenSum >= fiftypercent:
                    print_string += "N50: {}\n".format(int(val))
                    print_string += "L50: {}\n".format(count)
                    break
            print_string += "\n"

    print_string += lengthQual_table(df)

    if histogram:  # now make a histogram with read lengths
        histogram_lengths(df["length"], fastqBase.split('.')[0])
    print(print_string)


def lengthQual_table(df):
    """Create a table with lengths/basepairs on columns and qualities on rows"""
    # This block calculates the number of length bins for this data set
    lengthBinList = []
    size_map = [(1000, 250),
                (10000, 500),
                (40000, 1000),
                (100000, 5000),
                (500000, 20000),
                (1000000, 50000),
                (10000000000, 100000)]
    # first, figure out where we will start the table
    minlen = min(df["length"])
    current_val = 0
    firstDone = False
    for this_max_size, this_size_step in size_map:
        for this_bin in range(current_val, this_max_size, this_size_step):
            if minlen < this_bin:
                if not firstDone:
                    lengthBinList.append(prev)
                    firstDone = True
                lengthBinList.append(this_bin)
            prev = this_bin
        current_val = this_max_size
    # now figure out the largest bin
    maxLen = df["length"].max()
    first_index_gt_maxLen = next(i for i, v in enumerate(lengthBinList) if v > maxLen) + 1
    lengthBinList = lengthBinList[0:first_index_gt_maxLen]

    qualBinList = []
    increment_by = 1
    while len(qualBinList) == 0 or len(qualBinList) > 15:
        # now set up the bins for mean PHRED
        minQual = int(np.floor(min(df["meanQual"])))
        maxQual = int(np.ceil(max(df["meanQual"])))
        qualBinList = list(np.arange(minQual, maxQual + increment_by, increment_by))
        increment_by += 0.25

    # now make a table of read lengths
    bpTots = []
    readnumTots = []
    for row in range(len(lengthBinList)):
        dataNums = []
        readNums = []
        for column in range(len(qualBinList)):
            thisQuery = df.query("length >= {} and meanQual >= {}".format(
                lengthBinList[row], qualBinList[column]))
            dataNums.append(sum(thisQuery['length']))
            readNums.append(len(thisQuery.index))
        bpTots.append(dataNums)
        readnumTots.append(readNums)

    tables = {"Basepairs >= bin by mean PHRED and length": bpTots,
              "Number of reads >= bin by mean Phred+Len": readnumTots}
    print_table = ""
    for key in sorted(tables):
        # make a dataframe of our basepair distribution table
        dataDf = pd.DataFrame(tables[key], columns=["Q{}".format(x) for x in qualBinList])
        # add the min lengths as a column
        dataDf.insert(0, 'minLen', lengthBinList)
        print_table += pretty_print_table(dataDf, key)
    return print_table


def histogram_lengths(length, name_prefix):
    """Create a histogram of read counts per length."""
    counts = length.value_counts().to_frame(name="readCount")
    counts.index.rename('readLen', inplace=True)
    counts.sort_index(inplace=True)
    counts.to_csv("{}.hist.csv".format(name_prefix), index=True)


def pretty_print_table(df, title):
    print_string = ""
    dataframeStr = df.to_string(index=False)
    # this is the char width of the whole table printed
    lendataframeStr = len(dataframeStr.split('\n')[0])
    # this is the char width of the minlen portion of the printed table
    minLenLen = len(dataframeStr.split()[0])
    blank = " " * minLenLen
    # center the text on this offset as the table title
    txtoffset = lendataframeStr - minLenLen
    print_string += "\n{}{:^{offset}}\n".format(
        blank, title, offset=txtoffset)
    print_string += dataframeStr + "\n"
    return print_string


def run(args):
    """This just opens the fastq file and passes the info to the stats() function.
    This is a wrapper function that is accessed by pauvre_main.
    Useful since we can call stats() independently from other pauvre programs."""
    df = parse_fastq_length_meanqual(args.fastq)
    stats(df, args.fastq, histogram=args.histogram)
