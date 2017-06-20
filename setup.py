#!/usr/bin/env python3

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

# Tutorials on how to setup python package here:
#   - http://python-packaging.readthedocs.io/en/latest/testing.html
#   - https://jeffknupp.com/blog/2013/08/16/open-sourcing-a-python-project-the-right-way/

import os
from setuptools import setup, find_packages

version_py = os.path.join(os.path.dirname(__file__), 'pauvre', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','').strip()
print(version)


setup(name='pauvre',
      requires=['python (>3.0)'],
      version=version,
      description='Tools for plotting Oxford Nanopore and other long-read data.',
      long_description="""
          'pauvre' is a package for plotting Oxford Nanopore and other long read data.
          The name means 'poor' in French, a play on words to the oft-used 'pore' prefix
          for similar packages. The current implementation of this package only plots
          fastq quality. This package was designed for python 3, but it might work in
          python 2. You can visit the gitub page for more detailed information here:
          https://github.com/conchoecia/pauvre
      """,

      url='https://github.com/conchoecia/gloTK',
      author='Darrin Schultz',
      author_email='dts@ucsc.edu',
      classifiers=[
            'Development Status :: 2 - Pre-Alpha',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.5',
            'Operating System :: POSIX :: Linux',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Intended Audience :: Science/Research'
          ],
      license='GPLv3',
      provides=['pauvre'],
      packages=find_packages(),
      install_requires=[
          "matplotlib",
          "biopython",
          "pandas"
      ],
      entry_points = {
            'console_scripts': ['pauvre=pauvre.pauvre_main:main'],
        },
      zip_safe=False,
      include_package_data=True)
