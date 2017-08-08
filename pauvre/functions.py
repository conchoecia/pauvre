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

from Bio import SeqIO
import numpy as np

def parse_fastq_length_meanqual(fastq):
    """
    arguments:
     <fastq> the fastq file path. Hopefully it has been verified to exist already

    purpose:
     This function parses a fastq
    """
    length = []
    meanQual = []
    with open(fastq,"r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if len(record) > 0 and np.mean(record.letter_annotations["phred_quality"]) > 0:
                meanQual.append(_arithmetic_mean(record.letter_annotations["phred_quality"]))
                length.append(len(record))
    return length, meanQual

def _geometric_mean(phred_values):
    """in case I want geometric mean in the future, can calculate it like this"""
    #np.mean(record.letter_annotations["phred_quality"]))
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
    return np.power(10, (-1 * phred_values)/10)

def _erate_to_phred(erate_values):
    """
    converts a list or numpy array of error rates to a numpy array
    of phred values
    """
    if not isinstance(erate_values, np.ndarray):
        phred_values = np.array(erate_values)
    return -10 * np.log10(erate_values)
