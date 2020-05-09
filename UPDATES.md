## Updates:
- 20200508 - v0.2.2   - Minor changes to custommargin to make it work with newer matplotlib
- 20200501 - v0.2.1   - Minor changes related to setting matplotlib params
- 20200215 - v0.1.924 - Made some minor updates to work with python 3.7 and the latest version of pandas,
- 20171130 - v0.1.86 - some changes by @wdecoster to integrate `pauvre` into [nanoplot](https://github.com/wdecoster/NanoPlot),
  as well as some formatting changes that *may* make `pauvre` work better with python2.7. Adding Travis-CI functionality.
- 20171025 - v0.1.83 - added some changes to make marginplot interface
  with @wdecoster's [nanoPlot](https://github.com/wdecoster/NanoPlot)
  package, and made `pauvre stats` only output data tables for
  filtered reads. `pauvre stats` also now has the `--filt_maxlen`,
  `--filt_maxqual`, `--filt_minlen`, and `--filt_minqual` options.
- 20171018 - v0.1.8 - you can now filter reads and adjust the plotting viewing window.
  [See below for a demonstration.](#filter-reads-and-adjust-viewing-window) I added the following options:

```
  --filt_maxlen FILT_MAXLEN
                        This sets the max read length filter reads.
  --filt_maxqual FILT_MAXQUAL
                        This sets the max mean read quality to filter reads.
  --filt_minlen FILT_MINLEN
                        This sets the min read length to filter reads.
  --filt_minqual FILT_MINQUAL
                        This sets the min mean read quality to filter reads.
  --plot_maxlen PLOT_MAXLEN
                        Sets the maximum viewing area in the length dimension.
  --plot_maxqual PLOT_MAXQUAL
                        Sets the maximum viewing area in the quality
                        dimension.
  --plot_minlen PLOT_MINLEN
                        Sets the minimum viewing area in the length dimension.
  --plot_minqual PLOT_MINQUAL
                        Sets the minimum viewing area in the quality
                        dimension.
```
- 20171014 - uploading information on `pauvre redwood` and `pauvre synteny` usage.
- 20171012 - made `pauvre stats` more consistently produce useful histograms.
  `pauvre stats` now also calculates some statistics for different size ranges.
- 20170529 - added automatic scaling to the input fastq file. It
  scales to show the highest read quality and the top 99th percentile
  of reads by length.

