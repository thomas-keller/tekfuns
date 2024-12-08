
# tekfuns

<!-- badges: start -->
<!-- badges: end -->

The goal of tekfuns is to bundle together some clusterprofiler and deseq2 plots for rnaseq that I commonly need to make, hopefully with decent defaults. There are mostly sharp edges right now, so user beware.

## Installation

You can install the development version of tekfuns like so:

``` r
devtools::install('tekeller/tekfuns')
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(tekfuns)
# filter a gene set enrichment object by foldchange, return a cnetplot
# expects entrez ids
cf=cnetfilter(gse)

#difference fold change cutoff
cf=cnetfilter(gse,fcut=1.5)


#required input is a DESeqDataset object
#rnaplots will do the rest
#bundle of deseq and clusterprofiler plots
rp=rnaplots(dds,folder='exfolder',fprefix='exres')
#fcut and pcut are options now, for adjust cnetfilt and gsea, respectively
rp=rnaplots(dds,fcut=1.5,pcut=.1,folder='exfolder',fprefix='exres')
```

