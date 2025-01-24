# scRepertoire
<!-- badges: start -->
[![BioC status](http://www.bioconductor.org/shields/build/release/bioc/scRepertoire.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/scRepertoire)
[![R-CMD-check](https://github.com/ncborcherding/scRepertoire/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ncborcherding/scRepertoire/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/ncborcherding/scRepertoire/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ncborcherding/scRepertoire?branch=master)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://www.borch.dev/uploads/screpertoire/)
<!-- badges: end -->

## A toolkit for single-cell immune profiling

<img align="right" src="https://www.borch.dev/uploads/screpertoire/reference/figures/screpertoire_hex.png" width="305" height="352">

### Introduction
Single-cell sequencing is an emerging technology in the field of immunology and oncology that allows researchers to couple RNA quantification and other modalities, like immune cell receptor profiling at the level of an individual cell. A number of workflows and software packages have been created to process and analyze single-cell transcriptomic data. These packages allow users to take the vast dimensionality of the data generated in single-cell-based experiments and distill the data into novel insights. Unlike the transcriptomic field, there is a lack of options for software that allow for single-cell immune receptor profiling. Enabling users to easily combine RNA and immune profiling, the scRepertoire framework supports use of 10x, AIRR, BD, MiXCR, Omniscope, TRUST4, and WAT3R single-cell clonal formats and interaction with popular R-based single-cell data pipelines. 

### Applying Deep Learning to VDJ data
scRepertoire is compatible and integrated with the R packages [Trex](https://github.com/ncborcherding/Trex) for deep-learning-based autoencoding of the T cell receptor and [Ibex](https://github.com/ncborcherding/Ibex) for the B cell receptor. If you are interested in making your own deep-learning models with immune receptors, please see [immApex](https://github.com/ncborcherding/immApex).

## Working with scRepertoire

scRepertoire has a comprehensive [website](https://www.borch.dev/uploads/screpertoire/) for detailed tutorials and function information.

## Installation

### Installation of Master Branch

```R
devtools::install_github("ncborcherding/scRepertoire")
```

### Installing from Bioconductor
The current version of scRepertoire is also available on [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/scRepertoire.html).

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scRepertoire")
```

### Legacy Version 1
If you are looking for version 1 of scRepertoire - there is a static version available below:

```R
devtools::install_github("ncborcherding/scRepertoire@v1")
```

## Getting Data

Unfortunately, Github limits the size of individual files. In order to access the Seurat object paired with scRepertoire please download the .rda from [here](https://drive.google.com/file/d/1_YuRraDyg8UgF3oasjF0-jgPnwox-B24/view?usp=share_link).

## Bug Reports/New Features

#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/ncborcherding/scRepertoire/issues) with details of the issue.

- If possible please include a [reproducible example](https://reprex.tidyverse.org/). 
Alternatively, an example with the internal **scRep_example** and **contig_list** would 
be extremely helpful.

#### Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/ncborcherding/scRepertoire/issues).

#### [Pull Requests](https://github.com/ncborcherding/scRepertoire/pulls) are welcome for bug fixes, new features, or enhancements.

## Please Cite
scRepertoire 2 preprint is available at: [https://www.biorxiv.org/content/10.1101/2024.12.31.630854v1](https://www.biorxiv.org/content/10.1101/2024.12.31.630854v1)

*Yang, Q, Safina, K., Borcherding, N. (2024). scRepertoire 2: Enhanced and Efficient Toolkit for Single-Cell Immune Profiling https://doi.org/10.1101/2024.12.31.630854*

If you are building your own tool based on scRepertoire, reach out, we are happy to help and make things compatible.
