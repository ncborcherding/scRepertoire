#' Startrac-based diversity indices for single-cell RNA-seq
#'
#' @description This function utilizes the Startrac approach derived from 
#' [PMID: 30479382](https://pubmed.ncbi.nlm.nih.gov/30479382/).
#' Required to run the function, the "type" variable needs to include the 
#' difference in where the cells were derived. The output of this function 
#' will produce 3 indices: **expa** (clonal expansion), **migra** 
#' (cross-tissue migration), and **trans** (state transition). In order 
#' to understand the underlying analyses of the outputs please 
#' read and cite the linked manuscript. 
#' 
#' @examples
#' # Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' # Getting a sample of a Seurat object
#' scRep_example  <- get(data("scRep_example"))
#' scRep_example  <- combineExpression(combined, scRep_example)
#' scRep_example$Patient <- substring(scRep_example$orig.ident,1,3)
#' scRep_example$Type <- substring(scRep_example$orig.ident,4,4) 
#' 
#' # Using StartracDiversity()
#' StartracDiversity(scRep_example, 
#'                   type = "Type", 
#'                   group.by = "Patient")
#'
#' @param sc.data The single-cell object after [combineExpression()].
#' For SCE objects, the cluster variable must be in the meta data under 
#' "cluster".
#' @param cloneCall How to call the clone - VDJC gene (**gene**), 
#' CDR3 nucleotide (**nt**), CDR3 amino acid (**aa**),
#' VDJC gene + CDR3 nucleotide (**strict**) or a custom variable 
#' in the data. 
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL".
#' @param type The variable in the meta data that provides tissue type.
#' @param group.by The variable in the meta data to group by, often samples.
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any [hcl.pals][grDevices::hcl.pals].
#'
#' @export
#' @concept SC_Functions
#' @return ggplot object of Startrac diversity metrics
#' @author Liangtao Zheng
StartracDiversity <- function(sc.data,
                              cloneCall = "strict", 
                              chain = "both",
                              type = NULL,
                              group.by = NULL, 
                              exportTable = FALSE, 
                              palette = "inferno") {
    majorCluster <- NULL
    df <- .grabMeta(sc.data)
    cloneCall <- .theCall(df, cloneCall)
    barcodes <- rownames(df)
    colnames(df)[ncol(df)] <- "majorCluster"
    
    if (is.null(group.by)) {
       if(!"orig.ident" %in% colnames(df)) {
         stop("Please select a group.by variable")
       }
       group.by <- "orig.ident"
    }
    group.levels <- unique(df[,group.by])
    
    if (chain != "both") {
      df <- .offTheChain(df, chain, cloneCall)
    }

    df <- df %>%
              group_by(df[,group.by], df[,cloneCall]) %>%
              dplyr::mutate(n = n()) %>%
              as.data.frame()
                
    rownames(df) <- barcodes
    remove.pos <- which(df[,cloneCall] %in% c("", NA))
    if (length(remove.pos) > 0) {
      df <- df[-remove.pos,]
    }
    df[,"clone.status"] <- ifelse(df[,"n"] > 1, "Yes", "No")

    processed <- data.frame(rownames(df), df[,cloneCall], df$clone.status, 
                                df[,group.by], df[,"majorCluster"], 
                                df[,type], stringsAsFactors = FALSE)
    colnames(processed) <- c("Cell_Name", "clone.id", "clone.status", "patient", 
                                 "majorCluster", "loc")
    processed[processed == "NA"] <- NA
    processed <- na.omit(processed)
    
    mat.list <- list()
    for(i in seq_along(group.levels)) {
        mat.list[[i]] <- as.data.frame(.calIndex(processed[processed[,4] == group.levels[i],]))
    }
    names(mat.list) <- group.levels
    mat <- bind_rows(mat.list, .id = "group")
    if (exportTable) { 
      return(mat) 
    } 
    
    mat_melt <- melt(mat, id = c("group", "majorCluster"))
    values <-  .alphanumericalSort(unique(mat_melt$majorCluster))
    mat_melt$majorCluster <- factor(mat_melt$majorCluster, levels = values)
    mat_melt$value <- as.numeric(mat_melt$value)
    col <- length(unique(mat_melt[,"majorCluster"]))
        
    plot <- ggplot(mat_melt, aes(x=majorCluster, y=value)) +
                geom_boxplot(aes(fill = majorCluster), outlier.alpha = 0, na.rm = TRUE) +
                facet_grid(variable ~.) +
                theme_classic() +
                ylab("Index Score") +
                guides(fill="none") +
                scale_fill_manual(values = .colorizer(palette, col)) +
                theme(axis.title.x = element_blank())
      
    return(plot)
}


# Calculate cluster level indices
#' @importFrom plyr llply
.calIndex <- function(processed){
    clonotype.dist.cluster <- table(processed[,c("clone.id","majorCluster")])
    clonotype.dist.loc <- unclass(table(processed[,c("clone.id","loc")]))
    .entropy <- .mcol.entropy(clonotype.dist.cluster)
    .entropy.max <- log2(colSums(clonotype.dist.cluster > 0))
    expa <- 1-.entropy/.entropy.max

    .entropy.migr.max <- 1
    .entropy.tran.max <- 1
    
    clonotype.data <- 
      data.frame("clone.id"=rownames(clonotype.dist.loc),
                  "migr"=.mrow.entropy(clonotype.dist.loc)/.entropy.migr.max,
                  "tran"=.mrow.entropy(clonotype.dist.cluster)/.entropy.tran.max)
    weights.mtx <- sweep(clonotype.dist.cluster,2,colSums(clonotype.dist.cluster),"/")
    calIndex.matrix <- t(weights.mtx) %*% (as.matrix(clonotype.data[,c("migr","tran")]))
    calIndex.matrix <- cbind(calIndex.matrix, expa = expa)
    calIndex.matrix[is.nan(calIndex.matrix)] <- NA
    calIndex.matrix <- cbind(rownames(calIndex.matrix), calIndex.matrix)
    colnames(calIndex.matrix)[1] <- "majorCluster"
    return(calIndex.matrix)
}

# entropy of each row of the input matrix
.mrow.entropy <- function(x) {
    freqs <- sweep(x,1,rowSums(x),"/")
    H <- - rowSums(ifelse(freqs>0,freqs* log2(freqs),0))
    return(H)
}

# entropy of each column of the input matrix
.mcol.entropy <- function(x) {
    freqs <- sweep(x,2,colSums(x),"/")
    H <- - colSums(ifelse(freqs>0,freqs* log2(freqs),0))
    return(H)
}
