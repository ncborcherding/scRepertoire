# scRepertoire <img src="man/figures/screpertoire_hex.png" align="right" width="150"/>
<!-- badges: start -->
[![BioC status](http://www.bioconductor.org/shields/build/release/bioc/scRepertoire.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/scRepertoire)
[![R-CMD-check](https://github.com/ncborcherding/scRepertoire/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ncborcherding/scRepertoire/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/ncborcherding/scRepertoire/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ncborcherding/scRepertoire?branch=master)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://www.borch.dev/uploads/vignette/vignette)
<!-- badges: end -->

### Introduction
Single-cell sequencing is an emerging technology in the field of immunology and oncology that allows researchers to couple RNA quantification and other modalities, like immune cell receptor profiling at the level of an individual cell. A number of workflows and software packages have been created to process and analyze single-cell transcriptomic data. These packages allow users to take the vast dimensionality of the data generated in single-cell-based experiments and distill the data into novel insights. Unlike the transcriptomic field, there is a lack of options for software that allow for single-cell immune receptor profiling. Enabling users to easily combine RNA and immune profiling, the scRepertoire framework supports use of 10x, AIRR, BD, MiXCR, Omniscope, TRUST4, and WAT3R single-cell clonal formats and interaction with popular R-based single-cell data pipelines. 

### Applying Deep Learning to VDJ data
scRepertoire is compatible and integrated with the R packages [Trex](https://github.com/ncborcherding/Trex) for convolutional-neural-network-based autoencoding of the T cell receptor and [Ibex](https://github.com/ncborcherding/Ibex) for the B cell receptor. If you are interested in making your own deep-learning models with immune receptors, please see [immApex](https://github.com/ncborcherding/immApex).

## Installation

#### Installation of Master Branch
```
devtools::install_github("ncborcherding/scRepertoire")
```

### Installing from Bioconductor
The current version of scRepertoire is also available on [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/scRepertoire.html).

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scRepertoire")
```

#### Legacy Version 1
If you are looking for version 1 of scRepertoire - there is a static version available below:
```
devtools::install_github("ncborcherding/scRepertoire@v1")
```

## Getting Data

Unfortunately, Github limits the size of individual files. In order to access the full Seurat object paired with scRepertoire please download the .rda from [here](https://drive.google.com/file/d/1_YuRraDyg8UgF3oasjF0-jgPnwox-B24/view?usp=share_link). There is a smaller version of the cohort built into scRepertoire called **scRep_example**.

## Please Cite

We are working on submitting the scRepertoire v2 as a peer review article. For now, please cite scRepertoire using the doi from [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/scRepertoire.html).

*Borcherding, N, Yang, Q, & Safina, K. (2024). scRepertoire v2: A toolkit for single-cell immune receptor profiling. DOI: 10.18129/B9.bioc.scRepertoire*

## Bug Reports/New Features

#### [GitHub issue](https://github.com/ncborcherding/scRepertoire/issues) 

- If possible please include a [reproducible example](https://reprex.tidyverse.org/). 
Alternatively, an example with the internal **scRep_example** and **contig_list** would 
be extremely helpful.

#### [Features or Enhancements Requests](https://github.com/ncborcherding/scRepertoire/issues)

#### [Pull Requests](https://github.com/ncborcherding/scRepertoire/pulls) 
All are welcome for bug fixes, new features, or enhancements.

- Please set PR to merge with “dev” branch and provide description of what the PR 
contains (referencing existing issue(s) if appropriate).
  
## Contact
Questions, comments, suggestions, please feel free to reach out via the [GitHub Repository](https://github.com/ncborcherding/scRepertoire).

