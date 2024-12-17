
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

This is a basic example demonstrating a filtered cnet based on fold change, and and automatically generating some common deseq and clusterprofiler figures.
gsea is performed on hallmark and C2 pathways right now.

input: dds should be a DESEQdata object, with condition in the metadata (that is automatically the regression term)
given a prefix, it will make the following pdf figures:
cnetfilt- filtered cnetplot
enrichdot - gsea pathway dotplot
emap - gsea pathway term network plot
comp - gsea compare cluster (K=2)
gprof - gprofiler2 enrichment plot
pca -pca of genes colored by condition
heat
disp
topgdot

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

