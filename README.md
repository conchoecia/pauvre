## pauvre: a plotting package designed for nanopore reads

This package currently hosts one script for plotting.

- `pauvre marginplot`
  - takes a fastq file as input and outputs a marginal histogram with a heatmap.
  - transparency is the default

# Requirements

- You must have the following installed on your system to install this software:
  - python 3.x
  - matplotlib
  - pandas

# Installation

- Instructions to install on your mac or linux system. Not sure on
  Windows! Make sure python 3 is the active environment before
  installing.
  - `git clone https://github.com/conchoecia/pauvre.git`
  - `cd ./pauvre`
  - `pip install .`
- Or, install with pip
  - `pip install pauvre`

# Usage
- `marginplot`
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

![transparency](files/transparency.001.jpeg)

# Contributors

@conchoecia (Darrin Schultz)
