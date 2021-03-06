% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seuratFunctions.R
\name{combineExpression}
\alias{combineExpression}
\title{Adding clonotype information to a seurat or SCE object}
\usage{
combineExpression(
  df,
  sc,
  cloneCall = "gene+nt",
  groupBy = "none",
  proportion = TRUE,
  cloneTypes = c(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded
    = 1),
  filterNA = FALSE
)
}
\arguments{
\item{df}{The product of CombineTCR() or CombineBCR().}

\item{sc}{The seurat or SingleCellExperiment (SCE) object to attach}

\item{cloneCall}{How to call the clonotype - CDR3 gene (gene), 
CDR3 nucleotide (nt) CDR3 amino acid (aa), or 
CDR3 gene+nucleotide (gene+nt).}

\item{groupBy}{The column label in the combined contig object in which 
clonotype frequency will be calculated.}

\item{proportion}{Whether to use the total frequency (FALSE) or the 
proportion (TRUE) of the clonotype based on the groupBy variable.}

\item{cloneTypes}{The bins for the grouping based on frequency}

\item{filterNA}{Method to subset seurat object of barcodes without 
clonotype information}
}
\value{
seurat or SingleCellExperiment object with attached clonotype 
information
}
\description{
This function adds the immune receptor information to the seurat or 
SCE object to the meta data. By defualt this function also calculates 
the frequencies of the clonotypes by sequencing run (groupBy = "none"). 
To change how the frequencies are calculated, select a column header for 
the groupBy variable. Importantly, before using combineExpression() 
ensure the barcodes of the seurat or SCE object match the barcodes in the 
output of the combinedContig() call. Check changeNames() to change the 
prefix of the seurat object. If the dominant clonotypes have a greater 
frequency than 500, adjust the cloneTypes variable.
}
\examples{
#Getting the combined contigs
combined <- combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
rep(c("P", "T"), 3), cells ="T-AB")

#Getting a sample of a Seurat object
screp_example <- get(data("screp_example"))
sce <- suppressMessages(Seurat::UpdateSeuratObject(screp_example))
sce <- Seurat::as.SingleCellExperiment(sce)

#Using combineExpresion()
sce <- combineExpression(combined, sce)

}
