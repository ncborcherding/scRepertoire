# scRepertoire
#### A toolkit for single-cell immune profiling

<img align="right" src="https://github.com/ncborcherding/ncborcherding.github.io/blob/master/hex_sticker_v2.png" width="305" height="352">

### Introduction
Single-cell sequencing is an emerging technology in the field of immunology and oncology that allows researchers to couple RNA quantification and other modalities, like immune cell receptor profiling at the level of an individual cell. A number of workflows and software packages have been created to process and analyze single-cell transcriptomic data. These packages allow users to take the vast dimensionality of the data generated in single-cell-based experiments and distill the data into novel insights. Unlike the transcriptomic field, there is a lack of options for software that allow for single-cell immune receptor profiling. Enabling users to easily combine RNA and immune profiling, scRepertoire was built to process data derived from the 10x Genomics Chromium Immune Profiling for both T-cell receptor (TCR) and immunoglobulin (Ig) enrichment workflows and subsequently interacts with the popular Seurat R package. 

### R Packages Required
+  colorRamps
+  dplyr
+  ggalluvial
+  ggplots2
+  RColorBrewer
+  reshape2
+  Seurat
+  vegan

### Installation

```devtools::install_github("ncborcherding/scRepertoire")```

### Getting Data

Unfortunately, Github limits the size of individual files. In order to access the seurat object paired with scRepetoire please download the .rda from [here](https://drive.google.com/open?id=1np-EzG7U9W_Fz_SchBrsAhtqE3_rB_H9).

### Learning To Use scRepertoire

Vignette available [here](https://ncborcherding.github.io/vignettes/vignette.html), includes almost 13,000 T cells across 3 patients with renal clear cell carcinoma. 

### Contact
Questions, comments, suggestions, please feel free to contact Nick Borcherding via this repository, [email](mailto:ncborch@gmail.com), or using [twitter](https://twitter.com/theHumanBorch). 

