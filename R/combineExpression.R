#' Adding clonotype information to a single-cell object
#'
#' This function adds the immune receptor information to the Seurat or 
#' SCE object to the meta data. By default this function also calculates 
#' the frequencies and proportion of the clonotypes by sequencing 
#' run (**group.by** = "none"). To change how the frequencies/proportions
#' are calculated, select a column header for the **group.by** variable. 
#' Importantly, before using \code{\link{combineExpression}} ensure the 
#' barcodes of the single-cell object object match the barcodes in the output 
#' of the \code{\link{combineTCR}} or \code{\link{combineBCR}}. 
#'
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' #Getting a sample of a Seurat object
#' scRep_example <- get(data("scRep_example"))
#' sce <- suppressMessages(Seurat::UpdateSeuratObject(scRep_example))
#' sce <- Seurat::as.SingleCellExperiment(sce)
#' 
#' #Using combineExpresion()
#' sce <- combineExpression(combined, sce)
#' 
#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}} or a list of 
#' both c(\code{\link{combineTCR}}, c\code{\link{combineBCR}})
#' @param sc The seurat or SingleCellExperiment (SCE) object to attach
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt) CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by The column label in the combined clones in which 
#' clonotype frequency will be calculated. "none" will keep the list as is, 
#' while NULL will merge all the clones into a single data frame. 
#' @param proportion Whether to proportion (**TRUE**) or total frequency (**FALSE**) 
#' of the clonotype based on the group.by variable. 
#' @param cloneSize The bins for the grouping based on proportion or frequency. 
#' If proportion is **FALSE** and the *cloneSizes* are not set high enough based 
#' on frequency, the upper limit of *cloneSizes* will be automatically amended. 
#' @param filterNA Method to subset seurat object of barcodes without 
#' clonotype information
#' @param addLabel This will add a label to the frequency header, allowing
#' the user to try multiple group.by variables or recalculate frequencies after 
#' subsetting the data.
#' @importFrom dplyr bind_rows %>% summarise
#' @importFrom  rlang %||%
#' @importFrom SummarizedExperiment colData<- colData
#' @export
#' @return seurat or SingleCellExperiment object with attached clonotype 
#' information
#' 
combineExpression <- function(df, 
                              sc, 
                              cloneCall ="strict", 
                              chain = "both", 
                              group.by ="none", 
                              proportion = TRUE, 
                              filterNA = FALSE,
                              cloneSize = c(Rare = 1e-4,Small = 0.001,Medium = 0.01,Large = 0.1,Hyperexpanded = 1),
                              addLabel = FALSE) {
    call_time <- Sys.time()
  
    options( dplyr.summarise.inform = FALSE )
    if (!proportion & any(cloneSize < 1)) {
        stop("Adjust the cloneSize parameter - there are groupings < 1")
    }
    cloneSize <- c(None = 0, cloneSize)
    if (chain != "both") {
      df[[i]] <- .off.the.chain(df[[i]], chain, cloneCall)
    }
    df <- .checkList(df)
    cloneCall <- .theCall(cloneCall)
    
    #Getting Summaries of clones from combineTCR() or combineBCR()
    Con.df <- NULL
    meta <- .grabMeta(sc)
    cell.names <- rownames(meta)
    if (group.by == "none" | !is.null(group.by)) {
        for (i in seq_along(df)) {
      
            data <- data.frame(df[[i]], stringsAsFactors = FALSE)
            data2 <- unique(data[,c("barcode", cloneCall)])
            data2 <- na.omit(data2[data2[,"barcode"] %in% cell.names,])
            data2 <- data2 %>% 
                        group_by(data2[,cloneCall]) %>%
                        summarise(clonalProportion = n()/nrow(data2), 
                                  clonalFrequency = n())
            colnames(data2)[1] <- cloneCall
            data <- merge(data, data2, by = cloneCall, all = TRUE)
            data <- data[,c("barcode", "CTgene", "CTnt", 
                             "CTaa", "CTstrict", "clonalProportion", 
                             "clonalFrequency")]
            Con.df <- rbind.data.frame(Con.df, data)
        }
    } else if (group.by != "none" | !is.null(group.by)) {
        data <- data.frame(bind_rows(df), stringsAsFactors = FALSE)
        data2 <- na.omit(unique(data[,c("barcode", cloneCall, group.by)]))
        data2 <- data2[data2[,"barcode"] %in% cell.names, ]
        data2 <- as.data.frame(data2 %>% 
                                  group_by(data2[,cloneCall], data2[,group.by]) %>% 
                                  summarise(clonalProportion = n()/nrow(data2), 
                                            clonalFrequency = n())
        )
        
        colnames(data2)[c(1,2)] <- c(cloneCall, group.by)
        data <- merge(data, data2, by = cloneCall, all = TRUE)
        Con.df <- data[,c("barcode", "CTgene", "CTnt", 
                          "CTaa", "CTstrict", "clonalProportion", 
                          "clonalFrequency")]
        }
    #Detect if largest cloneSize category is too small for experiment and amend
    #this prevents a ton of NA values in the data
    if(!proportion && max(na.omit(Con.df[,"clonalFrequency"])) > cloneSize[length(cloneSize)]) {
      cloneSize[length(cloneSize)] <- max(na.omit(Con.df[,"clonalFrequency"]))
    }
    
    
    #Creating the bins for cloneSize
    Con.df$cloneSize <- NA
    for (x in seq_along(cloneSize)) { 
      names(cloneSize)[x] <- paste0(names(cloneSize[x]), ' (', cloneSize[x-1], 
        ' < X <= ', cloneSize[x], ')') 
    }
    
    if(proportion) {
      c.column <- "clonalProportion"
    } else {
      c.column <- "clonalFrequency"
    }
    #Assigning cloneSize
    for (i in 2:length(cloneSize)) { 
      Con.df$cloneSize <- ifelse(Con.df[,c.column] > cloneSize[i-1] & 
                                 Con.df[,c.column] <= cloneSize[i], 
                                 names(cloneSize[i]), 
                                 Con.df$cloneSize)
      }
    
    #Formating the meta data to add
    PreMeta <- unique(Con.df[,c("barcode", "CTgene", "CTnt", 
                "CTaa", "CTstrict", "clonalProportion", 
                "clonalFrequency", "cloneSize")])
    dup <- PreMeta$barcode[which(duplicated(PreMeta$barcode))]
    PreMeta <- PreMeta[PreMeta$barcode %!in% dup,]
    barcodes <- PreMeta$barcode
    PreMeta <- PreMeta[,-1]
    rownames(PreMeta) <- barcodes
    if (group.by != "none" && addLabel) {
      location <- which(colnames(PreMeta) %in% c("clonalProportion", 
                          "clonalFrequency"))
      colnames(PreMeta)[location] <- paste0(c("clonalProportion", 
                                            "clonalFrequency"), group.by)
    }
    

    
    if (is_seurat_object(sc)) { 
        if (length(which(rownames(PreMeta) %in% 
                         rownames(sc[[]])))/length(rownames(sc[[]])) < 0.01) {
          warning(.warn_str)
        }
        col.name <- names(PreMeta) %||% colnames(PreMeta)
        sc[[col.name]] <- PreMeta
    } else {
      rownames <- rownames(colData(sc))
      if (length(which(rownames(PreMeta) %in% 
                       rownames))/length(rownames) < 0.01) {
        warning(.warn_str) }
      colData(sc) <- cbind(colData(sc), PreMeta[rownames,])[, union(colnames(colData(sc)),  colnames(PreMeta))]
      rownames(colData(sc)) <- rownames  
    }
    if (filterNA) { sc <- .filteringNA(sc) }
    sc$cloneSize <- factor(sc$cloneSize, levels = rev(names(cloneSize)))
    
    if(is_seurat_object(sc)) {
        sc@commands[["combineExpression"]] <- make_screp_seurat_cmd(
            call_time, sc@active.assay
        )
    }
    sc
} 


.warn_str <- "< 1% of barcodes match: Ensure the barcodes in 
        the Seurat object match the barcodes in the combined immune receptor
        list from scRepertoire - most common issue is the addition of the 
        prefixes corresponding to 'samples' and 'ID' in the combineTCR/BCR() 
        functions"








