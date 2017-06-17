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


def stats(fastqName, length, meanQual):
    """
    arguments:
     <lengths> a list of read lengths
     <meanQual> a list of mean read qualities
     *NOTE* the indexing for each list is the same.
            element <lengths>[0] is the same read as <meanQual>[0]

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

    printString = "# Fastq stats for {}\n".format(fastqBase)

    printString += "numReads: {}\n".format(len(length))
    printString += "numBasepairs: {}\n".format(sum(length))
    printString += "meanLen: {}\n".format(np.mean(length))
    printString += "medianLen: {}\n".format(np.median(length))
    printString += "minLen: {}\n".format(min(length))
    maxLen = max(length)
    printString += "maxLen: {}\n".format(maxLen)

    #calculate the N50
    fiftypercent = 0.5 * sum(length)
    N50          = 0
    L50          = 0
    lenSum       = 0
    count        = 0
    for val in sorted(length, reverse=True):
        lenSum += val
        count += 1
        if lenSum >= fiftypercent:
            printString += "N50: {}\n".format(int(val))
            printString += "L50: {}\n".format(count)
            break

    #This block calculates the number of length bins for this data set
    scale = [10, 100, 1000, 10000, 100000, 1000000, 10000000]
    scaleChoice = 0
    for val in scale:
        thisDiv = maxLen / val
        if thisDiv < 10:
            scaleChoice = val
            break
    numBinSteps = int((int(maxLen / val) + 1) * 2)
    maxVal = int((int(maxLen / val) + 1) * val)
    binStepVal = int(maxVal/(numBinSteps))
    lengthBinList = []
    for i in range(numBinSteps):
        lengthBinList.append(binStepVal * i)


    #now set up the bins for mean PHRED
    maxQual = int(max(meanQual) + 1)
    stepInt = int(maxQual/5)
    rangeSet = set(range(0, maxQual+stepInt, stepInt))
    templist = [5, 10, 15, 20, 25, 30] + list(np.arange(29.5, 10, -4)) + list(np.arange(9.5, 0, -1))

    addlist = [x for x in templist if x < maxQual]
    for each in addlist:
        rangeSet.update([each])
        if len(rangeSet) == 10:
            break
    qualBinList = sorted(rangeSet)

    df = pd.DataFrame({'length': length,
                       'qual'  : meanQual})

    #now make a table of read lengths
    # row = j
    # column = i
    bpTots = []
    readnumTots = []
    for j in range(len(lengthBinList)):
        dataNums = []
        readNums = []
        for i in range(len(qualBinList)):
            thisQuery = df.query("length > {} and qual > {}".format(
                              lengthBinList[j], qualBinList[i]))
            dataNums.append(sum(thisQuery['length']))
            readNums.append(len(thisQuery.index))
        bpTots.append(dataNums)
        readnumTots.append(readNums)

    tables = {"Basepairs >= bin by mean PHRED and length": bpTots,
              "Number of reads >= bin by mean Phred+Len": readnumTots}
    for key in sorted(tables):
        #make a dataframe of our basepair distribution table
        dataDf = pd.DataFrame(tables[key], columns=["Q{}".format(x) for x in qualBinList])
        # add the min lengths as a column
        dataDf.insert(0, 'minLen', lengthBinList)
        dataframeStr = dataDf.to_string(index=False)
        # this is the char width of the whole table printed
        lendataframeStr = len(dataframeStr.split('\n')[0])
        # this is the char width of the minlen portion of the printed table
        minLenLen =  len(dataframeStr.split()[0])
        blank = " " * minLenLen
        # center the text on this offset as the table title
        txtoffset = lendataframeStr - minLenLen
        printString += "\n{}{:^{offset}}\n".format(
            blank, key, offset=txtoffset)
        printString += dataframeStr + "\n"


    print(printString)

def run(parser, args):
    """This just opens the fastq file and passes the info to the stats() function.
    This is a wrapper function that is accessed by pauvre_main.
    Useful since we can call stats() independently from other pauvre programs."""
    length, meanQual = parse_fastq_length_meanqual(args.fastq)
    stats(args.fastq, length, meanQual)
