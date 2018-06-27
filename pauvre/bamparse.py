#!/usr/bin/env python
# -*- coding: utf-8 -*-

# pauvre - just a pore plotting package
# Copyright (c) 2016-2018 Darrin T. Schultz. All rights reserved.
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
import pysam
import pandas as pd
import os

class BAMParse():
    """This class reads in a sam/bam file and constructs a pandas
    dataframe of all the relevant information for the reads to pass on
    and plot.
    """
    def __init__(self, filename, chrid = None, start = None,
                 stop = None, doubled = None):
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
        if chrid == None:
            self.chrid = self.sambam.references[0]
        else:
            self.chrid = chrid
        self.refindex = self.sambam.references.index(self.chrid)
        self.seqlength = self.sambam.lengths[self.refindex]
        self.true_seqlength = self.seqlength if not self.doubled else int(self.seqlength/2)
        if start == None or stop == None:
            self.start = 1
            self.stop = self.true_seqlength

        self.features = self.parse()
        self.features.sort_values(by=['POS','MAPLEN'], ascending=[True, False] ,inplace=True)
        self.features.reset_index(inplace=True)
        self.features.drop('index', 1, inplace=True)

        self.raw_depthmap = self.get_depthmap()
        self.features_depthmap = self.get_features_depthmap()

    def get_depthmap(self):
        depthmap = [0] * (self.stop - self.start + 1)
        for p in self.sambam.pileup(self.chrid, self.start, self.stop):
            index = p.reference_pos
            if index >= self.true_seqlength:
                index -= self.true_seqlength
            depthmap[index] += p.nsegments
        return depthmap

    def get_features_depthmap(self):
        """this method builds a more accurate pileup that is
        based on if there is actually a mapped base at any
        given position or not. better for long reads and RNA"""
        depthmap = [0] * (self.stop - self.start + 1)
        print("depthmap is: {} long".format(len(depthmap)))
        for index, row in self.features.iterrows():
            thisindex = row["POS"] - self.start
            for thistup in row["TUPS"]:
                b_type = thistup[1]
                b_len = thistup[0]
                if b_type == "M":
                    for j in range(b_len):
                        #this is necessary to reset the index if we wrap
                        # around to the beginning
                        if self.doubled and thisindex == len(depthmap):
                            thisindex = 0
                        depthmap[thisindex] += 1
                        thisindex += 1
                elif b_type in ["S", "H", "I"]:
                    pass
                elif b_type in ["D", "N"]:
                    thisindex += b_len
                    #this is necessary to reset the index if we wrap
                    # around to the beginning
                    if self.doubled and thisindex >= len(depthmap):
                        thisindex = thisindex - len(depthmap)

        return depthmap

    def parse(self):
        data = {'POS': [], 'MAPQ': [], 'TUPS': [] }
        for read in self.sambam.fetch(self.chrid, self.start, self.stop):
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
         an int

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
