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

# I took this code from https://github.com/arq5x/poretools/. Check it out. - DTS

import os.path
import sys
import argparse

#logger
import logging
logger = logging.getLogger('poretools')

# pauvre imports
import pauvre.version

#This class is used in argparse to expand the ~. This avoids errors caused on
# some systems.
class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                os.path.abspath(os.path.expanduser(values)))


def run_subtool(parser, args):
    if args.command == 'marginplot':
        import pauvre.marginplot as submodule

    # run the chosen submodule.
    submodule.run(parser, args)

class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                        action="store_true",
                        dest="quiet")

def main():

    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(prog='pauvre', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version", help="Installed pauvre version",
                        action="version",
                        version="%(prog)s " + str(pauvre.version.__version__))
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)

    #########################################
    # create the individual tool parsers
    #########################################



    ##########
    # mnplot
    ##########
    parser_mnplot = subparsers.add_parser('marginplot',
                                        help='plot a marginal histogram of a fastq file')
    parser_mnplot.add_argument('--fastq',
                               metavar='FASTQ',
                               action=FullPaths,
                               help='The input FASTQ file.')
    parser_mnplot.add_argument('--transparent',
                               metavar='TRANSPARENT',
                               default=True,
                               help="""Not the TV show. Set this to false if you
                               don't want the default output to have
                               transparency. Default is on.""")
    parser_mnplot.add_argument('--title',
                               metavar='TITLE',
                               default='Read length vs mean quality.',
                               help="""This sets the title for the whole plot.
                               Use --title "Crustacean's DNA read quality"
                               if you need single quote or apostrophe
                               inside title.""")
    parser_mnplot.add_argument('--maxlen',
                               metavar='MAXLEN',
                               type=int,
                               help="""This sets the max read length to plot.""")
    parser_mnplot.add_argument('--maxqual',
                               metavar='MAXQUAL',
                               type=int,
                               help="""This sets the max mean read quality
                               to plot.""")
    parser_mnplot.add_argument('--lengthbin',
                               metavar='LENGTHBIN',
                               type=int,
                               help="""This sets the bin size to use for length.""")
    parser_mnplot.add_argument('--qualbin',
                               metavar='QUALBIN',
                               type=float,
                               help="""This sets the bin size to use for quality""")
    parser_mnplot.add_argument('--fileform',
                               dest='fileform',
                               metavar='STRING',
                               choices=['png','pdf', 'eps', 'jpeg', 'jpg',
                                        'pdf', 'pgf', 'ps', 'raw', 'rgba',
                                        'svg', 'svgz', 'tif', 'tiff'],
                               default=['png'],
                               nargs='+',
                               help='Which output format would you like? Def.=png')
    parser_mnplot.add_argument('--dpi',
                               metavar='dpi',
                               default=600,
                               type=int,
                               help="""Change the dpi from the default 600
                               if you need it higher""")
    parser_mnplot.set_defaults(func=run_subtool)

    #######################################################
    # parse the args and call the selected function
    #######################################################
    args = parser.parse_args()

    # If there were no args, print the help function
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    print(args)

    if args.quiet:
        logger.setLevel(logging.ERROR)

    try:
        args.func(parser, args)
    except IOError as e:
         if e.errno != 32:  # ignore SIGPIPE
             raise

if __name__ == "__main__":
    main()
