#!/usr/bin/env python
# -*- coding: utf-8 -*-

# pauvre - just a pore plotting package
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

import os
import subprocess as sp
import unittest

class libSeq_test_case(unittest.TestCase):
    """Tests that different features of synplot work correctly
    """
    def setUp(self):
        # Here we must safely make the test directory if it does
        #  not exist
        scriptdir = os.path.dirname(os.path.realpath(__file__))
        testoutputdir = os.path.join(scriptdir, "testresults")
        if not os.path.exists(testoutputdir):
            os.makedirs(testoutputdir)
            # now we see if the specific subdirectory for the program
            #  we are testing exists
        self.thisoutdir = os.path.join(testoutputdir, "synplot")
        if not os.path.exists(self.thisoutdir):
            os.makedirs(self.thisoutdir)
            # now we change our working directory to the
            #  specific subdirectory for this file

        self.aln_dir = os.path.join(scriptdir, "testdata/alignments/")
        self.gff1 = os.path.join(scriptdir, "testdata/gff_files/Bf201706.gff")
        self.gff2 = os.path.join(scriptdir, "testdata/gff_files/JN392469.gff")
        self.gff3 = os.path.join(scriptdir, "testdata/gff_files/NC016117.gff")


    def test_normal_plotting_scenario(self):
        """This verifies that the LibSeq class is constructed with all of the
        parameters that are present in the meraculous config files"""
        os.chdir(self.thisoutdir)
        thiscommand =    """pauvre synplot --aln_dir {0} \
                         --fileform pdf \
                         --gff_paths {1} {2} {3} \
                         --center_on COX1 \
                         --gff_labels "B. forskalii" "P. bachei" "M. leidyi" """.format(
                             self.aln_dir, self.gff1,
                             self.gff2, self.gff3)

        data = sp.Popen(thiscommand, shell = True,
                        stdout=sp.PIPE)
        print("The data is:", data.communicate()[0].decode())
        print("The return code is: ", data.returncode)
        self.assertEqual(0, int(data.returncode))
