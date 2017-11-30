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

# I took this code from https://github.com/arq5x/poretools/. Check it out. - DTS

import sys
import os.path
import argparse

# logger
import logging
logger = logging.getLogger('poretools')

# pauvre imports
import pauvre.version
# This class is used in argparse to expand the ~. This avoids errors caused on
# some systems.


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                os.path.abspath(os.path.expanduser(values)))


class FullPathsList(argparse.Action):
    """Expand user- and relative-paths when a list of paths is passed to the
    program"""

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                [os.path.abspath(os.path.expanduser(value)) for value in values])


def run_subtool(parser, args):
    if args.command == 'marginplot':
        import pauvre.marginplot as submodule
    elif args.command == 'redwood':
        import pauvre.redwood as submodule
    elif args.command == 'stats':
        import pauvre.stats as submodule
    elif args.command == 'synplot':
        import pauvre.synplot as submodule
    # run the chosen submodule.
    submodule.run(args)


class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                          action="store_true",
                          dest="QUIET")

def main():

    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(
        prog='pauvre', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version", help="Installed pauvre version",
                        action="version",
                        version="%(prog)s " + str(pauvre.version.__version__))
    subparsers = parser.add_subparsers(
        title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)

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
    parser_mnplot.add_argument('--filt_maxlen',
                               type=int,
                               help="""This sets the max read length filter reads.""")
    parser_mnplot.add_argument('--filt_maxqual',
                               type=float,
                               help="""This sets the max mean read quality
                               to filter reads.""")
    parser_mnplot.add_argument('--filt_minlen',
                               type=int,
                               default=0,
                               help="""This sets the min read length to
                               filter reads.""")
    parser_mnplot.add_argument('--filt_minqual',
                               type=float,
                               default=0,
                               help="""This sets the min mean read quality
                               to filter reads.""")
    parser_mnplot.add_argument('--plot_maxlen',
                               type=int,
                               help="""Sets the maximum viewing area in the
                               length dimension.""")
    parser_mnplot.add_argument('--plot_maxqual',
                               type=float,
                               help="""Sets the maximum viewing area in the
                               quality dimension.""")
    parser_mnplot.add_argument('--plot_minlen',
                               type=int,
                               default=0,
                               help="""Sets the minimum viewing area in the
                               length dimension.""")
    parser_mnplot.add_argument('--plot_minqual',
                               type=float,
                               default=0,
                               help="""Sets the minimum viewing area in the
                               quality dimension.""")
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
                               choices=['png', 'pdf', 'eps', 'jpeg', 'jpg',
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
    # redwood
    #############
    parser_rwplot = subparsers.add_parser('redwood',
                                          help='make a redwood plot from a bam file')
    parser_rwplot.add_argument('-M', '--main_bam',
                               metavar='mainbam',
                               action=FullPaths,
                               help='The input filepath for the bam file to plot')
    parser_rwplot.add_argument('-R', '--rnaseq_bam',
                               metavar='rnabam',
                               action=FullPaths,
                               help='The input filepath for the rnaseq bam file to plot')
    parser_rwplot.add_argument('--gff',
                               metavar='gff',
                               action=FullPaths,
                               help="""The input filepath for the gff annotation
                               to plot""")
    parser_rwplot.add_argument('--fileform',
                               dest='fileform',
                               metavar='STRING',
                               choices=['png', 'pdf', 'eps', 'jpeg', 'jpg',
                                        'pdf', 'pgf', 'ps', 'raw', 'rgba',
                                        'svg', 'svgz', 'tif', 'tiff'],
                               default=['png'],
                               nargs='+',
                               help='Which output format would you like? Def.=png')
    parser_rwplot.add_argument('--sort',
                               dest='sort',
                               choices=['ALNLEN', 'TRULEN', 'MAPLEN', 'POS'],
                               default='ALNLEN',
                               help="""What value to use to sort the order in
                               which the reads are plotted?""")
    parser_rwplot.add_argument('--small_start',
                               dest='small_start',
                               choices=['inside', 'outside'],
                               default='inside',
                               help="""This determines where the shortest of the
                               filtered reads will appear on the redwood plot:
                               on the outside or on the inside? The default
                               option puts the longest reads on the outside and
                               the shortest reads on the inside.""")
    parser_rwplot.add_argument('--query',
                               dest='query',
                               default=['ALNLEN >= 10000', 'MAPLEN < reflength'],
                               nargs='+',
                               help='Query your reads to change plotting options')
    parser_rwplot.add_argument('-d', '--doubled',
                               dest='doubled',
                               choices=['main', 'rnaseq'],
                               default=[],
                               nargs='+',
                               help="""If your fasta file was doubled to map
                               circularly, use this flag or you may encounter
                               plotting errors. Accepts multiple arguments.
                               'main' is for the sam file passed with --sam,
                               'rnaseq' is for the sam file passed with --rnaseq""")
    parser_rwplot.add_argument('--dpi',
                               metavar='dpi',
                               default=600,
                               type=int,
                               help="""Change the dpi from the default 600
                               if you need it higher""")
    parser_rwplot.add_argument('-I', '--interlace',
                               action='store_true',
                               default=False,
                               help="""Interlace the reads so the pileup plot
                                       looks better""")
    parser_rwplot.add_argument('-i', '--invert',
                               action='store_true',
                               default=False,
                               help="""invert the image so that it looks better
                               on a dark backgroun. DOESN'T DO ANYTHING.""")
    parser_rwplot.add_argument('-L', '--log',
                               action='store_true',
                               default=False,
                               help="""Plot the RNAseq track with a log scale""")
    parser_rwplot.set_defaults(func=run_subtool)

    #############
    # stats
    #############
    parser_stats = subparsers.add_parser('stats',
                                         help='outputs stats from a fastq file')
    parser_stats.add_argument('-f', '--fastq',
                              metavar='FASTQ',
                              action=FullPaths,
                              help='The input FASTQ file.')
    parser_stats.add_argument('-H', '--histogram',
                              action='store_true',
                              help="""Make a histogram of the read lengths and
                               saves it to a new file""")
    parser_stats.add_argument('--filt_maxlen',
                              type=int,
                              help="""This sets the max read length filter reads.""")
    parser_stats.add_argument('--filt_maxqual',
                              type=float,
                              help="""This sets the max mean read quality
                               to filter reads.""")
    parser_stats.add_argument('--filt_minlen',
                              type=int,
                              default=0,
                              help="""This sets the min read length to
                               filter reads.""")
    parser_stats.add_argument('--filt_minqual',
                              type=float,
                              default=0,
                              help="""This sets the min mean read quality
                               to filter reads.""")

    parser_stats.set_defaults(func=run_subtool)

    #############
    # synplot
    #############
    parser_synplot = subparsers.add_parser('synplot',
                                           help="""make a synteny plot from a gff
                                        file, protein alignment, and partition
                                        file""")
    parser_synplot.add_argument('--gff_paths',
                                metavar='gff_paths',
                                action=FullPathsList,
                                nargs='+',
                                help="""The input filepath for the gff annotation
                                to plot""")
    parser_synplot.add_argument('--gff_labels',
                                metavar='gff_labels',
                                type=str,
                                nargs='+',
                                help="""In case the gff names and sequence names
                                don't match, change the labels that will appear
                                over the text.""")
    parser_synplot.add_argument('--dpi',
                                metavar='dpi',
                                default=600,
                                type=int,
                                help="""Change the dpi from the default 600
                                if you need it higher""")
    parser_synplot.add_argument('--optimum_order',
                                action='store_true',
                                help="""If selected, this doesn't plot the
                                optimum arrangement of things as they are input
                                into gff_paths. Instead, it uses the first gff
                                file as the top-most sequence in the plot, and
                                reorganizes the remaining gff files to minimize
                                the number of intersections.""")
    parser_synplot.add_argument('--aln_dir',
                                metavar='aln_dir',
                                action=FullPaths,
                                help="""The directory where all the fasta
                                alignments are contained.""")
    parser_synplot.add_argument('--stop_codons',
                                action='store_true',
                                default=True,
                                help="""Performs some internal corrections if
                                the gff annotation includes the stop
                                codons in the coding sequences.""")
    parser_synplot.add_argument('--center_on',
                                type=str,
                                default=None,
                                help="""centers the plot around the gene that
                                you pass as an argument""")
    parser_synplot.add_argument('--start_with_aligned_genes',
                                action='store_true',
                                default=False,
                                help="""Minimizes the number of intersections
                                but only selects combos where the first gene in
                                each sequence is aligned.""")
    parser_synplot.add_argument('--fileform',
                                dest='fileform',
                                metavar='STRING',
                                choices=['png', 'pdf', 'eps', 'jpeg', 'jpg',
                                         'pdf', 'pgf', 'ps', 'raw', 'rgba',
                                         'svg', 'svgz', 'tif', 'tiff'],
                                default=['png'],
                                nargs='+',
                                help='Which output format would you like? Def.=png')
    parser_synplot.add_argument('-o', '--output-base-name',
                                dest='BASENAME',
                                help='Specify a base name for the output file('
                                's). The input file base name is the '
                                'default.')
    parser_synplot.add_argument('-n', '--no-transparent',
                                dest='TRANSPARENT',
                                action='store_false',
                                help="""Not the TV show. Specify this option if
                               you don't want a transparent background. Default
                               is on.""")
    parser_synplot.add_argument('--sandwich',
                                action='store_true',
                                default=False,
                                help="""Put an additional copy of the first gff
                                file on the bottom of the plot for comparison.""")

    parser_synplot.set_defaults(func=run_subtool)

    #######################################################
    # parse the args and call the selected function
    #######################################################

    args = parser.parse_args()

    # If there were no args, print the help function
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # If there were no args, but someone selected a program,
    #  print the program's help.
    commandDict = {'redwood': parser_rwplot.print_help,
                   'marginplot': parser_mnplot.print_help,
                   'stats': parser_stats.print_help}

    if len(sys.argv) == 2:
        commandDict[args.command]()

        sys.exit(1)

    if args.QUIET:
        logger.setLevel(logging.ERROR)

    try:
        args.func(parser, args)
    except IOError as e:
        if e.errno != 32:  # ignore SIGPIPE
            raise


if __name__ == "__main__":
    main()
