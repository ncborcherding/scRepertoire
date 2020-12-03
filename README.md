# scRepertoire
#### A toolkit for single-cell immune profiling

<img align="right" src="https://github.com/ncborcherding/ncborcherding.github.io/blob/master/images/hex_sticker_v2.png" width="305" height="352">

### Introduction
Single-cell sequencing is an emerging technology in the field of immunology and oncology that allows researchers to couple RNA quantification and other modalities, like immune cell receptor profiling at the level of an individual cell. A number of workflows and software packages have been created to process and analyze single-cell transcriptomic data. These packages allow users to take the vast dimensionality of the data generated in single-cell-based experiments and distill the data into novel insights. Unlike the transcriptomic field, there is a lack of options for software that allow for single-cell immune receptor profiling. Enabling users to easily combine RNA and immune profiling, scRepertoire was built to process data derived from the 10x Genomics Chromium Immune Profiling for both T-cell receptor (TCR) and immunoglobulin (Ig) enrichment workflows and subsequently interacts with the popular Seurat R package. 

### Wrapper Functions
scRepertoire v1.0.2 has the functionality of the [powerTCR](https://github.com/hillarykoch/powerTCR) approach to comparing clone size distribution, [please cite](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006571) the manuscript if using the ```clonesizeDistribution()``` function. In addition, we recently added the [Startrac](https://github.com/Japrin/STARTRAC) clonotype metrics, if using the ```StartracDiversity()``` please read and cite [the accompanying article](https://www.nature.com/articles/s41586-018-0694-x).

### Installation

```
devtools::install_github("ncborcherding/scRepertoire")
```

**Note** Users with R versions < 4.0 or should use the dev version. Make sure you update tidyverse packages, like dplyr > 1.0. 

```
devtools::install_github("ncborcherding/scRepertoire@dev")
```

### Bioconductor

[scRepertoire](https://www.bioconductor.org/packages/release/bioc/html/scRepertoire.html) is now on Bioconductor. Use the link to find out more information. scRepertoire (for R v4 or greater) can be installed from Bioconductor with:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scRepertoire")
```

### Getting Data

Unfortunately, Github limits the size of individual files. In order to access the seurat object paired with scRepetoire please download the .rda from [here](https://drive.google.com/open?id=1np-EzG7U9W_Fz_SchBrsAhtqE3_rB_H9).

### Learning To Use scRepertoire

Vignette available [here](https://ncborcherding.github.io/vignettes/vignette.html), includes almost 13,000 T cells across 3 patients with renal clear cell carcinoma. 

### Please Cite
Feel free to use, edit, modify scRepertoire, but if you do, please cite the [manuscript](https://f1000research.com/articles/9-47/v1).

### Contact
Questions, comments, suggestions, please feel free to contact Nick Borcherding via this repository, [email](mailto:ncborch@gmail.com), or using [twitter](https://twitter.com/theHumanBorch). 

### More Information On VDJ Biology
Our publication is limited in scope and size and we were unable to review VDJ biology and issues relevant to high-throuhput immune receptor profiling. If you are interested in the biology of the process, please refer to the excellent reviews from [Jahong Ji](https://www.nature.com/articles/nri2941) and [Frederick Alt](https://www.cell.com/fulltext/S0092-8674(04)00039-X). Additionally, the complexities of sequencing and quantifications of VDJs are discussed [here.](https://www.nature.com/articles/nprot.2018.021?)

