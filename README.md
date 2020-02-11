# Tutorial on BLISS Downstream Analysis

Welcome to this tutorial. The aim is to provide a guide to the
preliminary downstream analysis of (s)BLISS data. A [Nature Protocols
article]() provides a summary of the sBLISS technology and a walkthrough
on the analysis of the results.

## Getting Started

### Requirements

In order to run this tutorial, the following tools are needed:
* Internet connection (to dynamically retrieve the ENSEMBL reference
annotation)
* R \>= 3.6.0
* Installed R packages: GenomicFeatures,
AnnotationHub, biomaRt, rtracklayer, data.table, circlize, ggplot2,
ggpubr, ggsci

### Setup

First, clone this repository on your computer and move to the new
directory:

``` bash
git clone git@github.com:BiCroLab/blissNPanalysis.git

cd blissNPanalysis
```

This tutorial assumes that the data required for its correct execution,
which can be downloaded from [here](), can be found in the *./data*
folder.

Once the data has been downloaded or moved to the data folder, open R
and continue to either of the following sections:

* [Human TK6 sBLISS analysis](README_human.md)
* [Mouse enterocyte sBLISS analysis](README_mouse.md)
