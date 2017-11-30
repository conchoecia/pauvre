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

from Bio import SeqIO
import copy
import numpy as np
import os
import pandas as pd
from sys import stderr
import os.path as opath
import matplotlib.pyplot as plt
import gzip


# this makes opening files more robust for different platforms
#  currently only used in GFFParse
import codecs

import warnings


def print_images(base_output_name, image_formats, dpi, path=None, transparent=False):
    file_base = opath.splitext(opath.basename(base_output_name))[0]
    for fmt in image_formats:
        if path:
            out_name = path
        else:
            out_name = "{}.{}".format(file_base, fmt)
        try:
            if fmt == 'png':
                plt.savefig(out_name, dpi=dpi, transparent=transparent)
            else:
                plt.savefig(out_name, format=fmt, transparent=transparent)
        except PermissionError:
            # thanks to https://github.com/wdecoster for the suggestion
            print("""You don't have permission to save pauvre plots to this
            directory. Try changing the directory and running the script again!""")


class GFFParse():
    def __init__(self, filename, stop_codons=None, species=None):
        self.filename = filename
        self.samplename = os.path.splitext(os.path.basename(filename))[0]
        self.species = species
        self.featureDict = {"name": [],
                            "featType": [],
                            "start": [],
                            "stop": [],
                            "strand": []}
        gffnames = ["sequence", "source", "featType", "start", "stop", "dunno1",
                    "strand", "dunno2", "tags"]
        self.features = pd.read_csv(self.filename, comment='#',
                                    sep='\t', names=gffnames)
        self.features['name'] = self.features['tags'].apply(self._get_name)
        self.features.drop('dunno1', 1, inplace=True)
        self.features.drop('dunno2', 1, inplace=True)
        self.features.reset_index(inplace=True, drop=True)
        # warn the user if there are CDS or gene entries not divisible by three
        self._check_triplets()
        # sort the database by start
        self.features.sort_values(by='start', ascending=True, inplace=True)
        if stop_codons:
            strip_codons = ['gene', 'CDS']
            # if the direction is forward, subtract three from the stop to bring it closer to the start
            self.features.loc[(self.features['featType'].isin(strip_codons)) & (self.features['strand'] == '+'), 'stop'] =\
                self.features.loc[(self.features['featType'].isin(strip_codons))
                                  & (self.features['strand'] == '+'), 'stop'] - 3
            # if the direction is reverse, add three to the start (since the coords are flip-flopped)
            self.features.loc[(self.features['featType'].isin(strip_codons)) & (self.features['strand'] == '-'), 'start'] =\
                self.features.loc[(self.features['featType'].isin(strip_codons))
                                  & (self.features['strand'] == '-'), 'start'] + 3
        self.features['center'] = self.features['start'] + \
            ((self.features['stop'] - self.features['start']) / 2)
        # we need to add one since it doesn't account for the last base otherwise
        self.features['width'] = abs(self.features['stop'] - self.features['start']) + 1
        self.features['lmost'] = self.features.apply(self._determine_lmost, axis=1)
        self.features['rmost'] = self.features.apply(self._determine_rmost, axis=1)
        self.features['track'] = 0
        if len(self.features.loc[self.features['tags'] == "Is_circular=true", 'stop']) < 1:
            raise IOError("""The GFF file needs to have a tag ending in "Is_circular=true"
            with a region from 1 to the number of bases in the mitogenome

            example:
            Bf201311	Geneious	region	1	13337	.	+	0	Is_circular=true
            """)
        self.seqlen = int(self.features.loc[self.features['tags'] == "Is_circular=true", 'stop'])
        self.features.reset_index(inplace=True, drop=True)
        #print("float", self.features.loc[self.features['name'] == 'COX1', 'center'])
        #print("float cat", len(self.features.loc[self.features['name'] == 'CAT', 'center']))
        # print(self.features)
        # print(self.seqlen)

    def set_features(self, new_features):
        """all this does is reset the features pandas dataframe"""
        self.features = new_features

    def get_unique_genes(self):
        """This returns a series of gene names"""
        plottable = self.features.query(
            "featType != 'tRNA' and featType != 'region' and featType != 'source'")
        return set(plottable['name'].unique())

    def shuffle(self):
        """
        this returns a list of all possible shuffles of features.
        A shuffle is when the frontmost bit of coding + noncoding DNA up
         until the next bit of coding DNA is removed and tagged on the
         end of the sequence. In this case this process is represented by
         shifting gff coordinates.
        """
        shuffles = []
        # get the index of the first element
        # get the index of the next thing
        # subtract the indices of everything, then reset the ones that are below
        # zero
        done = False
        shuffle_features = self.features[self.features['featType'].isin(
            ['gene', 'rRNA', 'CDS', 'tRNA'])].copy(deep=True)
        # we first add the shuffle features without reorganizing
        # print("shuffle\n",shuffle_features)
        add_first = copy.deepcopy(self)
        add_first.set_features(shuffle_features)
        shuffles.append(add_first)
        # first gene is changed with every iteration
        first_gene = list(shuffle_features['name'])[0]
        # absolute first is the first gene in the original gff file, used to determine if we are done in this while loop
        absolute_first = list(shuffle_features['name'])[0]
        while not done:
            # We need to prevent the case of shuffling in the middle of
            #  overlapped genes. Do this by ensuring that the the start of
            #  end of first gene is less than the start of the next gene.
            first_stop = int(shuffle_features.loc[shuffle_features['name'] == first_gene, 'stop'])
            next_gene = ""
            for next_index in range(1, len(shuffle_features)):
                # get the df of the next list, if len == 0, then it is a tRNA and we need to go to the next index
                next_gene_df = list(
                    shuffle_features[shuffle_features['featType'].isin(['gene', 'rRNA', 'CDS'])]['name'])
                if len(next_gene_df) != 0:
                    next_gene = next_gene_df[next_index]
                    next_start = int(
                        shuffle_features.loc[shuffle_features['name'] == next_gene, 'start'])
                    print("looking at {}, prev_stop is {}, start is {}".format(
                        next_gene, first_stop, next_start))
                    #print(shuffle_features[shuffle_features['featType'].isin(['gene', 'rRNA', 'CDS'])])
                    # if the gene we're looking at and the next one don't overlap, move on
                    if first_stop < next_start:
                        break
            print("next_gene before checking for first is {}".format(next_gene))
            if next_gene == absolute_first:
                done = True
                break
            # now we can reset the first gene for the next iteration
            first_gene = next_gene
            shuffle_features = shuffle_features.copy(deep=True)
            # figure out where the next start point is going to be
            next_start = int(shuffle_features.loc[shuffle_features['name'] == next_gene, 'start'])
            print('next gene: {}'.format(next_gene))
            shuffle_features['start'] = shuffle_features['start'] - next_start + 1
            shuffle_features['stop'] = shuffle_features['stop'] - next_start + 1
            shuffle_features['center'] = shuffle_features['center'] - next_start + 1
            # now correct the values that are less than 0
            shuffle_features.loc[shuffle_features['start'] < 1,
                                 'start'] = shuffle_features.loc[shuffle_features['start'] < 1, 'start'] + self.seqlen
            shuffle_features.loc[shuffle_features['stop'] < 1, 'stop'] = shuffle_features.loc[shuffle_features['stop']
                                                                                              < 1, 'start'] + shuffle_features.loc[shuffle_features['stop'] < 1, 'width']
            shuffle_features['center'] = shuffle_features['start'] + \
                ((shuffle_features['stop'] - shuffle_features['start']) / 2)
            shuffle_features['lmost'] = shuffle_features.apply(self._determine_lmost, axis=1)
            shuffle_features['rmost'] = shuffle_features.apply(self._determine_rmost, axis=1)
            shuffle_features.sort_values(by='start', ascending=True, inplace=True)
            shuffle_features.reset_index(inplace=True, drop=True)
            new_copy = copy.deepcopy(self)
            new_copy.set_features(shuffle_features)
            shuffles.append(new_copy)
        print("len shuffles: {}".format(len(shuffles)))
        return shuffles

    def couple(self, other_GFF, this_y=0, other_y=1):
        """
        Compares this set of features to another set and generates tuples of
        (x,y) coordinate pairs to input into lsi
        """
        other_features = other_GFF.features
        coordinates = []
        for thisname in self.features['name']:
            othermatch = other_features.loc[other_features['name'] == thisname, 'center']
            if len(othermatch) == 1:
                this_x = float(self.features.loc[self.features['name']
                                                 == thisname, 'center'])  # /self.seqlen
                other_x = float(othermatch)  # /other_GFF.seqlen
                # lsi can't handle vertical or horizontal lines, and we don't
                #  need them either for our comparison. Don't add if equivalent.
                if this_x != other_x:
                    these_coords = ((this_x, this_y), (other_x, other_y))
                    coordinates.append(these_coords)
        return coordinates

    def _check_triplets(self):
        """This method verifies that all entries of featType gene and CDS are
        divisible by three"""
        genesCDSs = self.features.query("featType == 'CDS' or featType == 'gene'")
        not_trips = genesCDSs.loc[((abs(genesCDSs['stop'] - genesCDSs['start']) + 1) % 3) > 0, ]
        if len(not_trips) > 0:
            print_string = ""
            print_string += "There are CDS and gene entries that are not divisible by three\n"
            print_string += str(not_trips)
            warnings.warn(print_string, SyntaxWarning)

    def _get_name(self, tag_value):
        """This extracts a name from a single row in 'tags' of the pandas
        dataframe
        """
        try:
            if ";" in tag_value:
                name = tag_value[5:].split(';')[0]
            else:
                name = tag_value[5:].split()[0]
        except:
            name = tag_value
            print("Couldn't correctly parse {}".format(
                tag_value))
        return name

    def _determine_lmost(self, row):
        """Booleans don't work well for pandas dataframes, so I need to use apply
        """
        if row['start'] < row['stop']:
            return row['start']
        else:
            return row['stop']

    def _determine_rmost(self, row):
        """Booleans don't work well for pandas dataframes, so I need to use apply
        """
        if row['start'] < row['stop']:
            return row['stop']
        else:
            return row['start']


def parse_fastq_length_meanqual(fastq):
    """
    arguments:
     <fastq> the fastq file path. Hopefully it has been verified to exist already

    purpose:
     This function parses a fastq and returns a pandas dataframe of read lengths
     and read meanQuals.
    """
    # First try to open the file with the gzip package. It will crash
    #  if the file is not gzipped, so this is an easy way to test if
    #  the fastq file is gzipped or not.
    try:
        handle = gzip.open(fastq, "rt")
        length, meanQual = _fastq_parse_helper(handle)
    except:
        handle = open(fastq, "r")
        length, meanQual = _fastq_parse_helper(handle)

    handle.close()
    df = pd.DataFrame(list(zip(length, meanQual)), columns=['length', 'meanQual'])
    return df


def filter_fastq_length_meanqual(df, min_len, max_len,
                                 min_mqual, max_mqual):
    querystring = "length >= {0} and meanQual >= {1}".format(min_len, min_mqual)
    if max_len != None:
        querystring += " and length <= {}".format(max_len)
    if max_mqual != None:
        querystring += " and meanQual <= {}".format(max_mqual)
    print("Keeping reads that satisfy: {}".format(querystring), file=stderr)
    filtdf = df.query(querystring)
    #filtdf["length"] = pd.to_numeric(filtdf["length"], errors='coerce')
    #filtdf["meanQual"] = pd.to_numeric(filtdf["meanQual"], errors='coerce')
    return filtdf


def _fastq_parse_helper(handle):
    length = []
    meanQual = []
    for record in SeqIO.parse(handle, "fastq"):
        if len(record) > 0:
            meanQual.append(_arithmetic_mean(record.letter_annotations["phred_quality"]))
            length.append(len(record))
    return length, meanQual


def _geometric_mean(phred_values):
    """in case I want geometric mean in the future, can calculate it like this"""
    # np.mean(record.letter_annotations["phred_quality"]))
    pass


def _arithmetic_mean(phred_values):
    """
    Convert Phred to 1-accuracy (error probabilities), calculate the arithmetic mean,
    log transform back to Phred.
    """
    if not isinstance(phred_values, np.ndarray):
        phred_values = np.array(phred_values)
    return _erate_to_phred(np.mean(_phred_to_erate(phred_values)))


def _phred_to_erate(phred_values):
    """
    converts a list or numpy array of phred values to a numpy array
    of error rates
    """
    if not isinstance(phred_values, np.ndarray):
        phred_values = np.array(phred_values)
    return np.power(10, (-1 * (phred_values / 10)))


def _erate_to_phred(erate_values):
    """
    converts a list or numpy array of error rates to a numpy array
    of phred values
    """
    if not isinstance(erate_values, np.ndarray):
        phred_values = np.array(erate_values)
    return -10 * np.log10(erate_values)
