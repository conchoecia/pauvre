### Note for pauvre

This version of bioawk contains function that quickly calculates the arithmetic mean of the quality scores. 

### Introduction

Bioawk is an extension to [Brian Kernighan's awk][1], adding the support of
several common biological data formats, including optionally gzip'ed BED, GFF,
SAM, VCF, FASTA/Q and TAB-delimited formats with column names. It also adds a
few built-in functions and an command line option to use TAB as the
input/output delimiter. When the new functionality is not used, bioawk is
intended to behave exactly the same as the original BWK awk.

The original awk requires a YACC-compatible parser generator (e.g. Byacc or
Bison). Bioawk further depends on [zlib][zlib] so as to work with gzip'd files.
