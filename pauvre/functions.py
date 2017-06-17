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
                meanQual.append(np.mean(record.letter_annotations["phred_quality"]))
                length.append(len(record))
    return length, meanQual
