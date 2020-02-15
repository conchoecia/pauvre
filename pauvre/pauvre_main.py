#!/usr/bin/env python
# -*- coding: utf-8 -*-

# pauvre - just a pore plotting package
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

# I modeled this code on https://github.com/arq5x/poretools/. Check it out. - DTS

import sys
import os.path
import argparse

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
    if args.command == 'browser':
        import pauvre.browser as submodule
    elif args.command == 'custommargin':
        import pauvre.custommargin as submodule
    elif args.command == 'marginplot':
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
    # browser
    #############
    parser_browser = subparsers.add_parser('browser',
                                          help="""an adaptable genome
                                          browser with various track types""")
    parser_browser.add_argument('-c', '--chromosomeid',
                                metavar = "Chr",
                                dest = 'CHR',
                                type = str,
                                help = """The fasta sequence to observe.
                                Use the header name of the fasta file
                                without the '>' character""")
    parser_browser.add_argument('--dpi',
                                metavar='dpi',
                                default=600,
                                type=int,
                                help="""Change the dpi from the default 600
                                if you need it higher""")
    parser_browser.add_argument('--fileform',
                                dest='fileform',
                                metavar='STRING',
                                choices=['png', 'pdf', 'eps', 'jpeg', 'jpg',
                                         'pdf', 'pgf', 'ps', 'raw', 'rgba',
                                         'svg', 'svgz', 'tif', 'tiff'],
                                default=['png'],
                                nargs='+',
                                help='Which output format would you like? Def.=png')
    parser_browser.add_argument("--no_timestamp",
                                  action = 'store_true',
                                  help="""Turn off time stamps in the filename
                                  output.""")
    parser_browser.add_argument('-o', '--output-base-name',
                               dest='BASENAME',
                               help="""Specify a base name for the output file(
                                    s). The input file base name is the
                                    default.""")
    parser_browser.add_argument('--path',
                               type=str,
                               help="""Set an explicit filepath for the output.
                               Only do this if you have selected one output type.""")
    parser_browser.add_argument('-p', '--plot_commands',
                                dest='CMD',
                                nargs = '+',
                                help="""Write strings here to select what
                                to plot. The format for each track is:
                                <type>:<path>:<style>:<options>

                                To plot the reference, the format is:
                                ref:<style>:<options>

                                Surround each track string with double
                                quotes and a space between subsequent strings.
                                "bam:mybam.bam:depth" "ref:colorful"
                                """)
    parser_browser.add_argument('--ratio',
                                nargs = '+',
                                type = float,
                                default=None,
                                help="""Enter the dimensions (arbitrary units)
                                to plot the figure. For example a figure that is
                                seven times wider than tall is:
                                --ratio 7 1""")
    parser_browser.add_argument('-r', '--reference',
                                metavar='REF',
                                dest = 'REF',
                                action=FullPaths,
                                help='The reference fasta file.')
    parser_browser.add_argument('--start',
                                metavar='START',
                                dest = 'START',
                                type = int,
                                help="""The start position to observe on the
                                fasta file. Uses 1-based indexing.""")
    parser_browser.add_argument('--stop',
                                metavar='STOP',
                                dest = 'STOP',
                                type = int,
                                help="""The stop position to observe on the
                                fasta file. Uses 1-based indexing.""")
    parser_browser.add_argument('-T', '--transparent',
                                action='store_false',
                                help="""Specify this option if you DON'T want a
                                transparent background. Default is on.""")
    parser_browser.set_defaults(func=run_subtool)

    ################
    # custommargin
    ################
    parser_custmar = subparsers.add_parser('custommargin',
                                          help='plot custom marginal histograms of tab-delimited files')
    parser_custmar.add_argument('--dpi',
                               metavar='dpi',
                               default=600,
                               type=int,
                               help="""Change the dpi from the default 600
                               if you need it higher""")
    parser_custmar.add_argument('--fileform',
                               dest='fileform',
                               metavar='STRING',
                               choices=['png', 'pdf', 'eps', 'jpeg', 'jpg',
                                        'pdf', 'pgf', 'ps', 'raw', 'rgba',
                                        'svg', 'svgz', 'tif', 'tiff'],
                               default=['png'],
                               nargs='+',
                               help='Which output format would you like? Def.=png')
    parser_custmar.add_argument('-i', '--input_file',
                                action=FullPaths,
                                help="""A tab-separated file with a header row
                                of column names.""")
    parser_custmar.add_argument('--xcol',
                                type=str,
                                help="""The column name of the data to plot on
                                the x-axis""")
    parser_custmar.add_argument('--ycol',
                                type=str,
                                help="""The column name of the data to plot on
                                the y-axis""")
    parser_custmar.add_argument('-n', '--no_transparent',
                               action='store_false',
                               help="""Specify this option if
                               you don't want a transparent background. Default
                               is on.""")
    parser_custmar.add_argument("--no_timestamp",
                                  action = 'store_true',
                                  help="""Turn off time stamps in the filename
                                  output.""")
    parser_custmar.add_argument('-o', '--output_base_name',
                                default = None,
                               help='Specify a base name for the output file('
                                    's). The input file base name is the '
                                    'default.')
    parser_custmar.add_argument('--plot_max_y',
                               type=int,
                               help="""Sets the maximum viewing area in the
                               length dimension.""")
    parser_custmar.add_argument('--plot_max_x',
                               type=float,
                               help="""Sets the maximum viewing area in the
                               quality dimension.""")
    parser_custmar.add_argument('--plot_min_y',
                               type=int,
                               default=0,
                               help="""Sets the minimum viewing area in the
                               length dimension. Default value = 0""")
    parser_custmar.add_argument('--plot_min_x',
                               type=float,
                               default=0,
                               help="""Sets the minimum viewing area in the
                               quality dimension. Default value = 0""")
    parser_custmar.add_argument('--square',
                               action='store_true',
                               help="""changes the hexmap into a square
                               map quantized by ints.""")
    parser_custmar.add_argument('-t', '--title',
                               type = str,
                               help="""This sets the title for the whole plot.
                               Use --title "Crustacean's DNA read quality"
                               if you need single quote or apostrophe
                               inside title.""")
    parser_custmar.add_argument('--ybin',
                               type=int,
                               help="""This sets the bin size to use for length.""")
    parser_custmar.add_argument('--xbin',
                               type=float,
                               help="""This sets the bin size to use for quality""")
    parser_custmar.add_argument('--add-yaxes',
                               dest='Y_AXES',
                               action='store_true',
                               help='Add Y-axes to both marginal histograms.')
    parser_custmar.set_defaults(func=run_subtool)


    #############
    # marginplot
    #############
    parser_mnplot = subparsers.add_parser('marginplot',
                                          help='plot a marginal histogram of a fastq file')
    parser_mnplot.add_argument('--dpi',
                               metavar='dpi',
                               default=600,
                               type=int,
                               help="""Change the dpi from the default 600
                               if you need it higher""")
    parser_mnplot.add_argument('-f', '--fastq',
                               metavar='FASTQ',
                               action=FullPaths,
                               help='The input FASTQ file.')
    parser_mnplot.add_argument('--fileform',
                               dest='fileform',
                               metavar='STRING',
                               choices=['png', 'pdf', 'eps', 'jpeg', 'jpg',
                                        'pdf', 'pgf', 'ps', 'raw', 'rgba',
                                        'svg', 'svgz', 'tif', 'tiff'],
                               default=['png'],
                               nargs='+',
                               help='Which output format would you like? Def.=png')
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
    parser_mnplot.add_argument('--kmerdf',
                               type = str,
                               default = None,
                               help = """Pass the filename of a data matrix if
                               you would like to plot read length
                               versus number of kmers in that read. The data matrix
                               is a tab-separated text file with columns
                               "id length numks and kmers", where:
                               <id> = read id
                               <length> = the length of the read
                               <numks> = the number of canonical kmers in the read
                               <kmers> = a list representation of kmers ie ['GAT', 'GTA']""")
    parser_mnplot.add_argument('-n', '--no_transparent',
                               dest='TRANSPARENT',
                               action='store_false',
                               help="""Specify this option if
                               you don't want a transparent background. Default
                               is on.""")
    parser_mnplot.add_argument("--no_timestamp",
                                  action = 'store_true',
                                  help="""Turn off time stamps in the filename
                                  output.""")
    parser_mnplot.add_argument('-o', '--output_base_name',
                               dest='BASENAME',
                               help='Specify a base name for the output file('
                                    's). The input file base name is the '
                                    'default.')
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
    parser_mnplot.add_argument('-t', '--title',
                               metavar='TITLE',
                               default='Read length vs mean quality',
                               help="""This sets the title for the whole plot.
                               Use --title "Crustacean's DNA read quality"
                               if you need single quote or apostrophe
                               inside title.""")
    parser_mnplot.add_argument('-y', '--add-yaxes',
                               dest='Y_AXES',
                               action='store_true',
                               help='Add Y-axes to both marginal histograms.')
    parser_mnplot.set_defaults(func=run_subtool)

    #############
    # redwood
    #############
    parser_redwood = subparsers.add_parser('redwood',
                                          help='make a redwood plot from a bam file')
    parser_redwood.add_argument('-d', '--doubled',
                               dest='doubled',
                               choices=['main', 'rnaseq'],
                               default=[],
                               nargs='+',
                               help="""If your fasta file was doubled to map
                               circularly, use this flag or you may encounter
                               plotting errors. Accepts multiple arguments.
                               'main' is for the sam file passed with --sam,
                               'rnaseq' is for the sam file passed with --rnaseq""")
    parser_redwood.add_argument('--dpi',
                               metavar='dpi',
                               default=600,
                               type=int,
                               help="""Change the dpi from the default 600
                               if you need it higher""")
    parser_redwood.add_argument('--fileform',
                               dest='fileform',
                               metavar='STRING',
                               choices=['png', 'pdf', 'eps', 'jpeg', 'jpg',
                                        'pdf', 'pgf', 'ps', 'raw', 'rgba',
                                        'svg', 'svgz', 'tif', 'tiff'],
                               default=['png'],
                               nargs='+',
                               help='Which output format would you like? Def.=png')
    parser_redwood.add_argument('--gff',
                               metavar='gff',
                               action=FullPaths,
                               help="""The input filepath for the gff annotation
                               to plot""")
    parser_redwood.add_argument('-I', '--interlace',
                               action='store_true',
                               default=False,
                               help="""Interlace the reads so the pileup plot
                                       looks better. Doesn't work currently""")
    parser_redwood.add_argument('-i', '--invert',
                               action='store_true',
                               default=False,
                               help="""invert the image so that it looks better
                               on a dark background. DOESN'T DO ANYTHING.""")
    parser_redwood.add_argument('-L', '--log',
                               action='store_true',
                               default=False,
                               help="""Plot the RNAseq track with a log scale""")
    parser_redwood.add_argument('-M', '--main_bam',
                               metavar='mainbam',
                               action=FullPaths,
                               help="""The input filepath for the bam file to
                               plot. Ideally was plotted with a fasta file that
                               is two copies of the mitochondrial genome
                               concatenated. This should be long reads (ONT, PB)
                               and will be displayed in the interior of the
                               redwood plot.""")
    parser_redwood.add_argument("--no_timestamp",
                                  action = 'store_true',
                                  help="""Turn off time stamps in the filename
                                  output.""")
    parser_redwood.add_argument('-o', '--output-base-name',
                               dest='BASENAME',
                               help='Specify a base name for the output file('
                                    's). The input file base name is the '
                                    'default.')
    parser_redwood.add_argument('--query',
                               dest='query',
                               default=['ALNLEN >= 10000', 'MAPLEN < reflength'],
                               nargs='+',
                               help='Query your reads to change plotting options')
    parser_redwood.add_argument('-R', '--rnaseq_bam',
                               metavar='rnabam',
                               action=FullPaths,
                               help='The input filepath for the rnaseq bam file to plot')
    parser_redwood.add_argument('--small_start',
                               dest='small_start',
                               choices=['inside', 'outside'],
                               default='inside',
                               help="""This determines where the shortest of the
                               filtered reads will appear on the redwood plot:
                               on the outside or on the inside? The default
                               option puts the longest reads on the outside and
                               the shortest reads on the inside.""")
    parser_redwood.add_argument('--sort',
                               dest='sort',
                               choices=['ALNLEN', 'TRULEN', 'MAPLEN', 'POS'],
                               default='ALNLEN',
                               help="""What value to use to sort the order in
                               which the reads are plotted?""")
    parser_redwood.add_argument('--ticks',
                               type = int,
                               nargs = '+',
                               default = [0, 10, 100, 1000],
                               help="""Specify control for the number of ticks.""")
    parser_redwood.add_argument('-T', '--transparent',
                                action='store_false',
                                help="""Specify this option if you DON'T want a
                                transparent background. Default is on.""")
    parser_redwood.set_defaults(func=run_subtool)

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
    parser_synplot.add_argument('--aln_dir',
                                metavar='aln_dir',
                                action=FullPaths,
                                help="""The directory where all the fasta
                                alignments are contained.""")
    parser_synplot.add_argument('--center_on',
                                type=str,
                                default=None,
                                help="""Centers the plot around the gene that
                                you pass as an argument. For example, if there
                                is a locus called 'COI' in the gff file and in
                                the alignments directory, center using:
                                --center_on COI""")
    parser_synplot.add_argument('--dpi',
                                metavar='dpi',
                                default=600,
                                type=int,
                                help="""Change the dpi from the default 600
                                if you need it higher""")
    parser_synplot.add_argument('--fileform',
                                dest='fileform',
                                choices=['png', 'pdf', 'eps', 'jpeg', 'jpg',
                                         'pdf', 'pgf', 'ps', 'raw', 'rgba',
                                         'svg', 'svgz', 'tif', 'tiff'],
                                default=['png'],
                                nargs='+',
                                help='Which output format would you like? Def.=png')
    parser_synplot.add_argument('--gff_paths',
                                metavar='gff_paths',
                                action=FullPathsList,
                                nargs='+',
                                help="""The input filepath. for the gff annotation
                                to plot. Individual filepaths separated by spaces.
                                For example, --gff_paths sp1.gff sp2.gff sp3.gff""")
    parser_synplot.add_argument('--gff_labels',
                                metavar='gff_labels',
                                type=str,
                                nargs='+',
                                help="""In case the gff names and sequence names
                                don't match, change the labels that will appear
                                over the text.""")
    parser_synplot.add_argument("--no_timestamp",
                                  action = 'store_true',
                                  help="""Turn off time stamps in the filename
                                  output.""")
    parser_synplot.add_argument('--optimum_order',
                                action='store_true',
                                help="""If selected, this doesn't plot the
                                optimum arrangement of things as they are input
                                into gff_paths. Instead, it uses the first gff
                                file as the top-most sequence in the plot, and
                                reorganizes the remaining gff files to minimize
                                the number of intersections.""")
    parser_synplot.add_argument('-o', '--output-basename',
                                dest='BASENAME',
                                help='Specify a base name for the output file('
                                's). The input file base name is the '
                                'default.')
    parser_synplot.add_argument('--ratio',
                                nargs = '+',
                                type = float,
                                default=None,
                                help="""Enter the dimensions (arbitrary units)
                                to plot the figure. For example a figure that is
                                seven times wider than tall is:
                                --ratio 7 1""")
    parser_synplot.add_argument('--sandwich',
                                action='store_true',
                                default=False,
                                help="""Put an additional copy of the first gff
                                file on the bottom of the plot for comparison.""")
    parser_synplot.add_argument('--start_with_aligned_genes',
                                action='store_true',
                                default=False,
                                help="""Minimizes the number of intersections
                                but only selects combos where the first gene in
                                each sequence is aligned.""")
    parser_synplot.add_argument('--stop_codons',
                                action='store_true',
                                default=True,
                                help="""Performs some internal corrections if
                                the gff annotation includes the stop
                                codons in the coding sequences.""")
    parser_synplot.add_argument('-T', '--transparent',
                                action='store_false',
                                help="""Specify this option if you DON'T want a
                                transparent background. Default is on.""")
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
    commandDict = {'browser': parser_browser.print_help,
                   'custommargin': parser_custmar.print_help,
                   'marginplot': parser_mnplot.print_help,
                   'redwood': parser_redwood.print_help,
                   'stats': parser_stats.print_help,
                   'synplot': parser_synplot.print_help}

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
