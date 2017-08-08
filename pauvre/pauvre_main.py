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

import sys
import os.path
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
    elif args.command == 'deathstar':
        import pauvre.deathstar as submodule
    elif args.command == 'stats':
        import pauvre.stats as submodule

    # run the chosen submodule.
    submodule.run(args)

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

    #############
    # marginplot
    #############
    parser_mnplot = subparsers.add_parser('marginplot',
                                        help='plot a marginal histogram of a fastq file')
    parser_mnplot.add_argument('-f', '--fastq',
                               metavar='FASTQ',
                               action=FullPaths,
                               help='The input FASTQ file.')
    parser_mnplot.add_argument('-n', '--no-transparent',
                               dest='TRANSPARENT',
                               action='store_false',
                               help="""Not the TV show. Specify this option if
                               you don't want a transparent background. Default
                               is on.""")
    parser_mnplot.add_argument('-t', '--title',
                               metavar='TITLE',
                               default='Read length vs mean quality',
                               help="""This sets the title for the whole plot.
                               Use --title "Crustacean's DNA read quality"
                               if you need single quote or apostrophe
                               inside title.""")
    parser_mnplot.add_argument('--maxlen',
                               metavar='MAXLEN',
                               type=int,
                               help="""This sets the max read length to plot.""")
    parser_mnplot.add_argument('-m', '--maxqual',
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
    parser_mnplot.add_argument('-y', '--add-yaxes',
                               dest='Y_AXES',
                               action='store_true',
                               help='Add Y-axes to both marginal histograms.')
    parser_mnplot.add_argument('--fileform',
                               dest='fileform',
                               metavar='STRING',
                               choices=['png','pdf', 'eps', 'jpeg', 'jpg',
                                        'pdf', 'pgf', 'ps', 'raw', 'rgba',
                                        'svg', 'svgz', 'tif', 'tiff'],
                               default=['png'],
                               nargs='+',
                               help='Which output format would you like? Def.=png')
    parser_mnplot.add_argument('-o', '--output-base-name',
                               dest='BASENAME',
                               help='Specify a base name for the output file('
                                    's). The input file base name is the '
                                    'default.')
    parser_mnplot.add_argument('-d', '--dpi',
                               metavar='dpi',
                               default=600,
                               type=int,
                               help="""Change the dpi from the default 600
                               if you need it higher""")
    parser_mnplot.set_defaults(func=run_subtool)

    #############
    # deathstar
    #############
    parser_dsplot = subparsers.add_parser('deathstar',
                                        help='make a deathstar plot from a bam file')
    parser_dsplot.add_argument('-M','--main_bam',
                               metavar='mainbam',
                               action=FullPaths,
                               help='The input filepath for the bam file to plot')
    parser_dsplot.add_argument('-R','--rnaseq_bam',
                               metavar='rnabam',
                               action=FullPaths,
                               help='The input filepath for the rnaseq bam file to plot')
    parser_dsplot.add_argument('--gff',
                               metavar='gff',
                               action=FullPaths,
                               help="""The input filepath for the gff annotation
                               to plot""")
    parser_dsplot.add_argument('--fileform',
                               dest='fileform',
                               metavar='STRING',
                               choices=['png','pdf', 'eps', 'jpeg', 'jpg',
                                        'pdf', 'pgf', 'ps', 'raw', 'rgba',
                                        'svg', 'svgz', 'tif', 'tiff'],
                               default=['png'],
                               nargs='+',
                               help='Which output format would you like? Def.=png')
    parser_dsplot.add_argument('--sort',
                               dest='sort',
                               choices=['ALNLEN', 'TRULEN', 'MAPLEN', 'POS'],
                               default='ALNLEN',
                               help="""What value to use to sort the order in
                               which the reads are plotted?""")
    parser_dsplot.add_argument('--small_start',
                               dest='small_start',
                               choices=['inside', 'outside'],
                               default='inside',
                               help="""This determines where the shortest of the
                               filtered reads will appear on the deathstar plot:
                               on the outside or on the inside?""")
    parser_dsplot.add_argument('--query',
                               dest='query',
                               default=['ALNLEN >= 10000','MAPLEN < reflength'],
                               nargs='+',
                               help='Query your reads to change plotting options')
    parser_dsplot.add_argument('-d', '--doubled',
                               dest='doubled',
                               choices=['main','rnaseq'],
                               default=[],
                               nargs='+',
                               help="""If your fasta file was doubled to map
                               circularly, use this flag or you may encounter
                               plotting errors. Accepts multiple arguments.
                               'main' is for the sam file passed with --sam,
                               'rnaseq' is for the sam file passed with --rnaseq""")
    parser_dsplot.add_argument('--dpi',
                               metavar='dpi',
                               default=600,
                               type=int,
                               help="""Change the dpi from the default 600
                               if you need it higher""")
    parser_dsplot.add_argument('-I', '--interlace',
                               action='store_true',
                               default=False,
                               help="""Interlace the reads so the pileup plot
                                       looks better""")
    parser_dsplot.add_argument('-L', '--log',
                               action='store_true',
                               default=False,
                               help="""Plot the RNAseq track with a log scale""")
    parser_dsplot.set_defaults(func=run_subtool)

    #############
    # stats
    #############
    parser_stats = subparsers.add_parser('stats',
                                        help='outputs stats from a fastq file')
    parser_stats.add_argument('-f', '--fastq',
                               metavar='FASTQ',
                               action=FullPaths,
                               help='The input FASTQ file.')
    parser_stats.set_defaults(func=run_subtool)

    #######################################################
    # parse the args and call the selected function
    #######################################################

    args = parser.parse_args()

    # If there were no args, print the help function
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # If there were no args, but someone selected a program,
    #  print the program's help.
    commandDict={'deathstar' : parser_dsplot.print_help,
                 'marginplot': parser_mnplot.print_help,
                 'stats'     : parser_stats.print_help}

    if len(sys.argv)==2:
        commandDict[args.command]()

        sys.exit(1)

    if args.quiet:
        logger.setLevel(logging.ERROR)

    try:
        args.func(parser, args)
    except IOError as e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

if __name__ == "__main__":
    main()
