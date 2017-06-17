## pauvre: a plotting package designed for nanopore reads

This package currently hosts one script for plotting.

- `pauvre marginplot`
  - takes a fastq file as input and outputs a marginal histogram with a heatmap.
- transparency is the default

## Updates:
- 20170529 - added automatic scaling to the input fastq file. It
  scales to show the highest read quality and the top 99th percentile
  of reads by length.

# Requirements

- You must have the following installed on your system to install this software:
  - python 3.x
  - matplotlib
  - biopython
  - pandas

# Installation

- Instructions to install on your mac or linux system. Not sure on
  Windows! Make sure *python 3* is the active environment before
  installing.
  - `git clone https://github.com/conchoecia/pauvre.git`
  - `cd ./pauvre`
  - `pip3 install .`
- Or, install with pip
  - `pip3 install pauvre`

# Usage
- `marginplot`
  - Make the default plot showing the 99th percentile of longest reads
    - `pauvre marginplot --fastq miniDSMN15.fastq`
    - ![default](files/default_miniDSMN15.png)
  - Make a marginal histogram for ONT 2D or 1D^2 cDNA data with a
    lower maxlen and higher maxqual.
    - `pauvre marginplot --maxlen 4000 --maxqual 25 --lengthbin 50 --fileform pdf png --qualbin 0.5 --fastq miniDSMN15.fastq`
    - ![example1](files/miniDSMN15.png)
  - Plot ONT 1D data with a large tail
    - `pauvre marginplot --maxlen 100000 --maxqual 15 --lengthbin 500  <myfile>.fastq`
  - Get more resolution on lengths
    - `pauvre marginplot --maxlen 100000 --lengthbin 5  <myfile>.fastq`
  - Turn off transparency if you just want a white background
    - `pauvre marginplot --transparent False <myfile>.fastq`
    - Note: transparency is the default behavior
      - ![transparency](files/transparency.001.jpeg)

# Contributors

@conchoecia (Darrin Schultz)
