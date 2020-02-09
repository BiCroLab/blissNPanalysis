# Tutorial on BLISS Downstream Analysis

Welcome to this tutorial. The aim is to provide a guide to the
preliminary downstream analysis of BLISS data. A [Nature Protocols
article]() provides a summary of the BLISS technology and a walkthrough
on the analysis and interpretation of the results.

## Getting Started

### Requirements

In order to run this tutorial, the following tools are needed: \*
Internet connection (to dynamically retrieve the ENSEMBL reference
annotation) \* R \>= 3.6 \* Installed R packages: GenomicFeatures,
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

* [Human TK6 BLISS analysis](README_human.md)
* [Mouse CD73 BLISS analysis](README_mouse.md)
