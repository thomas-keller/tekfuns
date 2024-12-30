
# tekfuns

<!-- badges: start -->
<!-- badges: end -->

The goal of tekfuns is to bundle together some clusterprofiler and deseq2 plots for rnaseq that I commonly need to make, hopefully with decent defaults. There are mostly sharp edges right now, so user beware.

Notably, I added in some transcript-level analysis with fishpond and network/regulon prediction with genie3. If you enable these options, expect the analysis to take several more minutes. 

Right now with all options takes ~7-10 minutes.

also, to save individual output (all results are returned in a list object), provide a folder AND a file prefix.

running genie3 on the whole rnaseq dataset produces a very large file, shortening to top 200 genes for now. hardcoded for now, allow option later.

clusterprofiler has changed some of the option calls in the latest bioconductor. If you get an error, or things arent displaying correctly, update bioconductor and the associated libraries.

12/30/24

genie3 now predicts based on a list of [transcription factors](https://inesdesantiago.github.io/SeqQC.blog/TFlists/Final_TFlist.txt) . Thanks for the list and the [extensive post on how it was created](https://seqqc.wordpress.com/2020/12/05/where-to-find-a-comprehensive-list-of-potential-human-transcription-factors/)

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

#all bells and whistles, run transcript-level analysis and regulon inference
#se is a tximeta object from salmon
rp=rnaplots(dds,sw=se,regulons=TRUE,pcut=.1,nfcut=1,folder='exfolder',fprefix='exres')

```

