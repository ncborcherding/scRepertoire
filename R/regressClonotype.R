#' Normalize feature space of single-cell expression by clonotype
#' 
#' This function reduces the influence of clonotype on the feature space of the 
#' single-cell object. It first removes the genes that comprise the adaptive immune 
#' receptors from the list of variable genes. Next, it uses the 
#' \href{https://github.com/immunogenomics/harmony}{harmony R package} to re-integrate 
#' data to minimize effects of clonotype and other variables indicated by groupVariable.
#'
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' #Getting a sample of a Seurat object
#' screp_example <- get(data("screp_example"))
#' sce <- suppressMessages(Seurat::UpdateSeuratObject(screp_example))
#' 
#' #Using combineExpresion()
#' sce <- combineExpression(combined, sce)
#' 
#' sce[["RNA"]] <- sce[["integrated"]]
#' sce <- regressClonotype(sce, groupVariable = NULL, reduction = "UMAP", dims = 1:5)
#' 
#' @param sc The seurat or SCE object to visualize after combineExpression(). 
#' For SCE objects, the cluster variable must be in the meta data under 
#' "cluster".
#' @param cloneCall How to call the clonotype - CDR3 gene (gene), 
#' CDR3 nucleotide (nt) or CDR3 amino acid (aa), or 
#' CDR3 gene+nucleotide (gene+nt).
#' @param groupVariable Additional features to normalize the single cell object with
#' @param reduction Re-perform the dimensional reduction - "UMAP" or "TSNE"
#' @param dims The number of dimensions to use in the reduction (i.e., 1:20)
#' @param ... Other parameters to be called from Harmony package
#' @importFrom Seurat RunUMAP RunTSNE
#' @importFrom scater calculatePCA runUMAP runTSNE
#' @importFrom harmony RunHarmony HarmonyMatrix
#' @importFrom SingleCellExperiment reducedDim<-
#' @export
#' @return Seurat or SingleCellExperiment object 

regressClonotype <- function(sc, cloneCall="gene+nt", 
                             groupVariable = NULL, 
                             reduction = NULL, 
                             dims = NULL, 
                             ...) {
    cloneCall <- theCall(cloneCall)
    
    groupVariable <- c(groupVariable, cloneCall)
    unwanted_genes <- "^IGHV*|^IGHJ*|^IGHD*|^IGKV*|^IGLV*|^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*"
    if (inherits(x=sc, what ="Seurat")) {
        unwanted_genes <- grep(pattern = unwanted_genes, x = sc[["RNA"]]@var.features, value = TRUE)
        remove_genes <- sc[["RNA"]]@var.features %in% unwanted_genes
        sc[["RNA"]]@var.features = sc[["RNA"]]@var.features[!remove_genes]
        
        sc <- RunHarmony(sc, groupVariable, verbose = FALSE, ...)
        
        if(reduction == "UMAP") {
            sc <- RunUMAP(sc, reduction = "harmony", dims = dims)
        }
        if(reduction == "TSNE") {
            sc <- RunTSNE(sc, reduction = "harmony", dims = dims)
        }
    } else {
        meta_data <- grabMeta(sc)
        PCA <- calculatePCA(sc, subset_row = rownames(sc)[!grepl(unwanted_genes, rownames(sc))])
        
        my_harmony_embeddings <- HarmonyMatrix(
            PCA, meta_data, groupVariable,
            do_pca = FALSE, 
            verbose = FALSE, ...)
    reducedDim(sc) <- my_harmony_embeddings
    if(reduction == "UMAP") {
        sc <- runUMAP(sc, dimred="PCA", n_dimred=dims[length(dims)])
    }
    if(reduction == "TSNE") {
        sc <- runTSNE(sc, dimred="PCA", n_dimred=dims[length(dims)])
    }
    }
    return(sc)
}
