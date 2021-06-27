#' Adding clonotype information to a seurat or SCE object
#'
#' This function adds the immune receptor information to the seurat or 
#' SCE object to the meta data. By defualt this function also calculates 
#' the frequencies of the clonotypes by sequencing run (groupBy = "none"). 
#' To change how the frequencies are calculated, select a column header for 
#' the groupBy variable. Importantly, before using combineExpression() 
#' ensure the barcodes of the seurat or SCE object match the barcodes in the 
#' output of the combinedContig() call. Check changeNames() to change the 
#' prefix of the seurat object. If combining more 
#' than one immune receptor type, barcodes with both receptors will be removed
#' during the combination process. 
#'
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' #Getting a sample of a Seurat object
#' screp_example <- get(data("screp_example"))
#' sce <- suppressMessages(Seurat::UpdateSeuratObject(screp_example))
#' sce <- Seurat::as.SingleCellExperiment(sce)
#' 
#' #Using combineExpresion()
#' sce <- combineExpression(combined, sce)
#' 
#' @param df The product of CombineTCR() or CombineBCR() or a list of 
#' both c(CombineTCR(), combineBCR())
#' @param sc The seurat or SingleCellExperiment (SCE) object to attach
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt) CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (gene+nt).
#' @param groupBy The column label in the combined contig object in which 
#' clonotype frequency will be calculated.
#' @param proportion Whether to use the total frequency (FALSE) or the 
#' proportion (TRUE) of the clonotype based on the groupBy variable.
#' @param cloneTypes The bins for the grouping based on frequency
#' @param filterNA Method to subset seurat object of barcodes without 
#' clonotype information
#' @importFrom dplyr bind_rows %>% summarise
#' @importFrom  rlang %||%
#' @importFrom SummarizedExperiment colData<- colData
#' @export
#' @return seurat or SingleCellExperiment object with attached clonotype 
#' information
#' 

combineExpression <- function(df, sc, cloneCall="gene+nt", groupBy="none", 
                              proportion = TRUE,
                            cloneTypes=c(Rare = 1e-4, Small = 0.001, 
                            Medium = 0.01, Large = 0.1, Hyperexpanded = 1), filterNA = FALSE) {
  options( dplyr.summarise.inform = FALSE )
    cloneTypes <- c(None = 0, cloneTypes)
    df <- checkList(df)
    cloneCall <- theCall(cloneCall)
    Con.df <- NULL
    meta <- grabMeta(sc)
    cell.names <- rownames(meta)
    if (groupBy == "none") {
        for (i in seq_along(df)) {
            data <- data.frame(df[[i]], stringsAsFactors = FALSE)
            data2 <- unique(data[,c("barcode", cloneCall)])
            data2 <- na.omit(data2[data2[,"barcode"] %in% cell.names,])
            if (proportion == TRUE) {
                data2 <- data2 %>% group_by(data2[,cloneCall]) %>%
                    summarise(Frequency = n()/nrow(data2))
            } else {
            data2 <- data2 %>% group_by(data2[,cloneCall]) %>%
                summarise(Frequency = n())
            }
            colnames(data2)[1] <- cloneCall
            data <- merge(data, data2, by = cloneCall, all = TRUE)
            data <- data[,c("barcode", "CTgene", "CTnt", 
                             "CTaa", "CTstrict", "Frequency")]
            Con.df <- rbind.data.frame(Con.df, data)
        }
    } else if (groupBy != "none") {
        data <- data.frame(bind_rows(df), stringsAsFactors = FALSE)
        data2 <- na.omit(unique(data[,c("barcode", cloneCall, groupBy)]))
        data2 <- data2[data2[,"barcode"] %in% cell.names, ]
        data2 <- as.data.frame(data2 %>% group_by(data2[,cloneCall], 
                    data2[,groupBy]) %>% summarise(Frequency = n()))
        colnames(data2)[c(1,2)] <- c(cloneCall, groupBy)
        x <- unique(data[,groupBy])
        for (i in seq_along(x)) {
            sub1 <- subset(data, data[,groupBy] == x[i])
            sub2 <- subset(data2, data2[,groupBy] == x[i])
            merge <- merge(sub1, sub2, by=cloneCall)
            if (proportion == TRUE) {
                merge$Frequency <- merge$Frequency/length(merge$Frequency)
            }
            Con.df <- rbind.data.frame(Con.df, merge)
        } 
        nsize <- Con.df %>% group_by(Con.df[,paste0(groupBy, ".x")])  %>% summarise(n = n())
    }
    
    Con.df$cloneType <- NA
    for (x in seq_along(cloneTypes)) { names(cloneTypes)[x] <- 
        paste0(names(cloneTypes[x]), ' (', cloneTypes[x-1], 
        ' < X <= ', cloneTypes[x], ')') }
    for (i in 2:length(cloneTypes)) { Con.df$cloneType <- 
        ifelse(Con.df$Frequency > cloneTypes[i-1] & Con.df$Frequency 
        <= cloneTypes[i], names(cloneTypes[i]), Con.df$cloneType) }
    PreMeta <- unique(Con.df[,c("barcode", "CTgene", "CTnt", 
                "CTaa", "CTstrict", "Frequency", "cloneType")])
    `%!in%` = Negate(`%in%`)
    PreMeta <- PreMeta[PreMeta$barcode %!in% dup,]
    rownames(PreMeta) <- PreMeta$barcode
    if (inherits(x=sc, what ="Seurat")) { 
        if (length(which(rownames(PreMeta) %in% 
                         rownames(sc[[]])))/length(rownames(sc[[]])) < 0.01) {
          warning("< 1% of barcodes match: Ensure the barcodes in 
            the Seurat object match the 
            barcodes in the combined immune receptor list from 
            scRepertoire - most common issue is the addition of the 
            prefixes corresponding to `samples` and 'ID' in the combineTCR/BCR() 
            functions")
        }
        col.name <- names(PreMeta) %||% colnames(PreMeta)
        sc[[col.name]] <- PreMeta
    } else {
      rownames <- rownames(colData(sc))
      if (length(which(rownames(PreMeta) %in% 
                       rownames))/length(rownames) < 0.01) {
        warning("< 1% of barcodes match: Ensure the barcodes 
          in the SingleCellExperiment object match the 
          barcodes in the combined immune receptor list from 
          scRepertoire - most common issue is the addition of the 
          prefixes corresponding to `samples` and 'ID' in the combineTCR/BCR() 
          functions")
      colData(sc) <- cbind(colData(sc), PreMeta[rownames,])[, union(colnames(colData(sc)),  colnames(PreMeta))]
      rownames(colData(sc)) <- rownames  } 
    }
    if (filterNA == TRUE) { sc <- filteringNA(sc) }
    return(sc) 
}

#' Highlighting specific clonotypes in Seurat
#'
#' Use a specific clonotype sequence to highlight on top of the dimensional 
#' reduction in seurat object.
#'
#' @examples
#' #' #Getting the combined contigs
#' combined <- combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' #Getting a sample of a Seurat object
#' screp_example <- get(data("screp_example"))
#' 
#' #Using combineExpresion()
#' screp_example <- combineExpression(combined, screp_example )
#' 
#' #Using highlightClonotype()
#' screp_example  <- highlightClonotypes(screp_example, cloneCall= "aa", 
#' sequence = c("CAVNGGSQGNLIF_CSAEREDTDTQYF"))
#' 
#' @param sc The seurat object to attach
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (gene+nt).
#' @param sequence The specific sequence or sequence to highlight
#'
#' @export
#' @return DimPlot with highlighted clonotypes
highlightClonotypes <- function(sc, 
                            cloneCall = c("gene", "nt", "aa", "gene+nt"), 
                            sequence = NULL){
    if (!inherits(x=sc, what ="Seurat")) {
        stop("Object indicated is not of class 'Seurat', make sure you 
            are using the correct data.") }
    cloneCall <- theCall(cloneCall)
    meta <- grabMeta(sc)
    meta$highlight <- NA
    for(i in seq_along(sequence)) {
        meta$highlight <- ifelse(meta[,cloneCall] == sequence[i], 
                            paste("Clonotype", i, sep=""),  meta$highlight) }
    meta <- meta[,-(which(colnames(meta) == "cluster"))]
    col.name <- names(meta) %||% colnames(meta)
    sc[[col.name]] <- meta
    return(sc)
}

#' Exploring interaction of clonotypes by seurat or SCE dynamics
#'
#' View the proportional contribution of clonotypes by seurat or SCE object 
#' meta data after combineExpression(). The visualization is based on the 
#' ggalluvial package, which requires the aesthetics to be part of the axes 
#' that are visualized. Therefore, alpha, facet, and color should be part of 
#' the the axes you wish to view or will add an additional stratum/column to 
#' the end of the graph.
#'
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' #Getting a sample of a Seurat object
#' screp_example <- get(data("screp_example"))
#' sce <- suppressMessages(Seurat::UpdateSeuratObject(screp_example))
#' sce <- Seurat::as.SingleCellExperiment(sce)
#' 
#' #Using combineExpresion()
#' sce <- combineExpression(combined, sce)
#' 
#' #Using alluvialClonotypes()
#' alluvialClonotypes(sce, cloneCall = "gene", 
#' y.axes = c("Patient", "cluster"), color = "cluster")
#' 
#' @param sc The seurat or SCE object to visualize after combineExpression(). 
#' For SCE objects, the cluster variable must be in the meta data under 
#' "cluster".
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt) or CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (gene+nt).
#' @param y.axes The columns that will separate the proportional 
#' visualizations.
#' @param color The column header or clonotype(s) to be highlighted.
#' @param facet The column label to separate.
#' @param alpha The column header to have gradated opacity.
#' @import ggplot2
#' 
#' @importFrom ggalluvial StatStratum geom_flow geom_stratum to_lodes_form geom_alluvium
#' @importFrom dplyr %>% mutate
#' 
#' @export
#' @return Alluvial ggplot comparing clonotype distribution across 
#' selected parameters.
alluvialClonotypes <- function(sc, 
                                cloneCall = c("gene", "nt", "aa", "gene+nt"), 
                                y.axes = NULL, color = NULL, alpha = NULL, 
                                facet = NULL) {
    checkSingleObject(sc)
    cloneCall <- theCall(cloneCall)
    if (length(y.axes) == 0) {
        stop("Make sure you have selected the variable(s) to visualize") }
    meta <- grabMeta(sc)
    meta$barcodes <- rownames(meta)
    check <- colnames(meta) == color
    if (length(unique(check)) == 1 & unique(check)[1] == FALSE & 
        !is.null(color)) {
        meta <- meta %>% mutate("clonotype(s)" = ifelse(meta[,cloneCall] %in% 
            color, "Selected", "Other"))
        color <- "clonotype(s)" }
    
    y.axes <- unique(c(y.axes, color, alpha, facet))
    set.axes <- seq_along(y.axes)
    meta2 <- meta[,c(y.axes, color, alpha, facet, cloneCall, "barcodes")]
    meta2 <- unique(na.omit(meta2[!duplicated(as.list(meta2))]))

    lodes <- makingLodes(meta2, color, alpha, facet, set.axes) 
    plot <- ggplot(data = lodes, aes(x = x, stratum = stratum, 
                alluvium = alluvium, label = stratum)) +
        geom_stratum() + theme_classic() +
        theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())  
    if (is.null(color) & is.null(alpha)) {
        plot <- plot + geom_alluvium()
    } else if (!is.null(color) & is.null(alpha)) {
        plot <- plot+geom_flow(aes(fill = lodes[,color]), 
            stat = "alluvium", lode.guidance = "forward") + labs(fill = color)
    } else if (is.null(color) & !is.null(alpha)) {
        plot <- plot + geom_flow(aes(alpha = lodes[,alpha]), stat = "alluvium",
                lode.guidance = "forward") + labs(alpha = alpha)
    }else {
        plot <- plot+geom_flow(aes(alpha=lodes[,alpha], fill=lodes[,color]),
                        stat = "alluvium", lode.guidance = "forward") + 
                        labs(fill = color, alpha = alpha) }
    if (length(facet) == 1 & length(facet) < 2) {
        plot <- plot + facet_wrap(.~lodes[,facet], scales="free_y")
    } else if (length(facet) == 0) { plot <- plot }
    plot <- plot + geom_text(stat = ggalluvial::StatStratum, infer.label = FALSE, reverse = TRUE, size = 2)
    return(plot)}


#' Visualize the number of single cells with clonotype frequencies by cluster
#'
#' View the count of clonotypes frequency group in seurat or SCE object 
#' meta data after combineExpression(). The visualization will take the 
#' new meta data variable "cloneType" and plot the number of cells with
#' each designation using a secondary variable, like cluster. Credit to 
#' the idea goes to Drs. Carmona and Andreatta and their work with
#' \href{https://github.com/carmonalab/ProjecTILs}{ProjectTIL}.
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
#' #Using occupiedscRepertoire()
#' occupiedscRepertoire(sce, x.axis = "cluster")
#' table <- occupiedscRepertoire(sce, x.axis = "cluster", exportTable = TRUE)
#' 
#' @param sc The seurat or SCE object to visualize after combineExpression(). 
#' For SCE objects, the cluster variable must be in the meta data under 
#' "cluster".
#' @param x.axis The variable in the meta data to graph along the x.axis
#' @param label Include the number of clonotype in each category by x.axis variable
#' @param proportion Convert the stacked bars into relative proption
#' @param exportTable Exports a table of the data into the global 
#' environment in addition to the visualization
#' @importFrom dplyr %>% group_by mutate
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
#' @return Stacked bar plot of counts of cells by clonotype frequency group

occupiedscRepertoire <- function(sc, x.axis = "cluster", 
                                 label = TRUE, proportion = FALSE, 
                                 exportTable = FALSE) {
    checkSingleObject(sc)
    meta <- grabMeta(sc)
    meta <- melt(table(meta[!is.na(meta$Frequency), 
                c(x.axis, "cloneType")]), varnames = c(x.axis, "cloneType"))
    meta <- meta[meta$value != 0,]
    if(proportion == TRUE) {
      meta <- meta %>%
        group_by(meta[,1]) %>%
        mutate(total = sum(value), 
               prop = value/total)
      meta <- as.data.frame(meta)
    }
    if (exportTable == TRUE) {
        return(meta)
    }
    col <- length(unique(meta$cloneType))
    if(proportion == TRUE) {
      plot <- ggplot(meta, aes(x = meta[,x.axis], y = prop, fill = cloneType)) + 
        geom_bar(stat = "identity") 
         
    } else {
      plot <- ggplot(meta, aes(x = meta[,x.axis], y = value, fill = cloneType)) + 
        geom_bar(stat = "identity") 
      
    } 
    plot <- plot + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(values = c(colorblind_vector(col))) + 
      ylab("Single Cells") + 
      theme_classic() + 
      theme(axis.title.x = element_blank())
    if (label == TRUE) {
        plot <- plot + geom_text(aes(label = value), position = position_stack(vjust = 0.5))
      }
    return(plot)
}

#' Visualize distribution of clonal frequency overlaid on dimensional reduction plots
#'
#' This function allows the user to visualize the clonal expansion by overlaying the 
#' cells with specific clonal frequency onto the dimensional reduction plots in Seurat.
#' Credit to the idea goes to Dr. Carmona and his work with
#' \href{https://github.com/carmonalab/ProjecTILs}{ProjectTIL}.
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
#' #Using clonalOverlay()
#' clonalOverlay(sce, freq.cutpoint = 0.3, bins = 5) 
#' 
#' @param sc The seurat or SCE object to visualize after combineExpression(). 
#' @param reduction The dimensional reduction to visualize
#' @param freq.cutpoint The overlay cutpoint to include, this corresponds to the 
#' Frequency variable in the single-cell objecter
#' @param bins The number of contours to the overlay
#' @param facet meta data variable to facet the comparison
#' 
#' @import ggplot2
#' @importFrom SeuratObject Embeddings
#' @export
#' @author Francesco Mazziotta, Nick Borcherding
#' 
#' @return ggplot object

clonalOverlay <- function(sc, reduction = NULL, freq.cutpoint = 30, bins = 25, facet = NULL) {
  checkSingleObject(sc)
  if (is.null(reduction)) {
    tmp <- data.frame(sc[[]], identity = sc@active.ident, Embeddings(sc, reduction = "pca"))
  } else {
    tmp <- data.frame(sc[[]], identity = sc@active.ident, Embeddings(sc, reduction = reduction))
  }
  if (!is.null(facet)) {
    facet <- tmp[,facet]
    tmp <-data.frame(facet, tmp)
  }
  tmp$include <- ifelse(tmp$Frequency >= freq.cutpoint, "Yes", NA)
  tmp2 <- subset(tmp, include == "Yes")
  plot <- ggplot(tmp2, mapping = aes(x=tmp2[,(ncol(tmp2)-2)], y = tmp2[,(ncol(tmp2)-1)])) +
    geom_point(tmp, mapping = aes(x=as.numeric(tmp[,(ncol(tmp)-2)]), y = as.numeric(tmp[,(ncol(tmp)-1)]), color = tmp[,"identity"]), size= 0.5) +
    geom_density_2d(color = "black", lwd=0.25, bins = bins) + 
    theme_classic() +
    labs(color = "Active Identity") +
    xlab("Dimension 1") + ylab("Dimension 2")
  if (!is.null(facet)) {
    plot <- plot + facet_wrap(~facet) 
  }
  return(plot)
}
