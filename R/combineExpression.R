#' Adding clonotype information to a Seurat or SCE object
#'
#' This function adds the immune receptor information to the Seurat or 
#' SCE object to the meta data. By default this function also calculates 
#' the frequencies of the clonotypes by sequencing run (group.by = "none"). 
#' To change how the frequencies are calculated, select a column header for 
#' the group.by variable. Importantly, before using combineExpression() 
#' ensure the barcodes of the seurat or SCE object match the barcodes in the 
#' output of the combinedContig() call. Check changeNames() to change the 
#' prefix of the Seurat object. If combining more than one immune receptor 
#' type, barcodes with both receptors will be removed during the combination 
#' process.
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
#' @param df The product of CombineTCR() or CombineBCR() or a list of 
#' both c(CombineTCR(), combineBCR())
#' @param sc The seurat or SingleCellExperiment (SCE) object to attach
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt) CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by The column label in the combined contig object in which 
#' clonotype frequency will be calculated. "none" will keep the list as is, 
#' while NULL will merge all the contigs into a single data frame. 
#' @param proportion Whether to use the total frequency (FALSE) or the 
#' proportion (TRUE) of the clonotype based on the group.by variable.
#' @param cloneTypes The bins for the grouping based on frequency
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
combineExpression <- function(
    df, 
    sc, 
    cloneCall ="strict", 
    chain = "both", 
    group.by ="none", 
    proportion = TRUE, 
    filterNA = FALSE,
    cloneTypes = c(
        Rare = 1e-4,Small = 0.001,Medium = 0.01,Large = 0.1,Hyperexpanded = 1
    ),
    addLabel = FALSE
) {
    call_time <- Sys.time()
  
    options( dplyr.summarise.inform = FALSE )
    cloneTypes <- c(None = 0, cloneTypes)
    df <- checkList(df)
    cloneCall <- theCall(cloneCall)
    Con.df <- NULL
    meta <- grabMeta(sc)
    cell.names <- rownames(meta)
    if (group.by == "none" | !is.null(group.by)) {
        for (i in seq_along(df)) {
            if (chain != "both") {
                df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
            }
            data <- data.frame(df[[i]], stringsAsFactors = FALSE)
            data2 <- unique(data[,c("barcode", cloneCall)])
            data2 <- na.omit(data2[data2[,"barcode"] %in% cell.names,])
            if (proportion) {
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
    } else if (group.by != "none" | is.null(group.by)) {
        data <- data.frame(bind_rows(df), stringsAsFactors = FALSE)
        data2 <- na.omit(unique(data[,c("barcode", cloneCall, group.by)]))
        data2 <- data2[data2[,"barcode"] %in% cell.names, ]
        data2 <- as.data.frame(data2 %>% group_by(data2[,cloneCall], 
                    data2[,group.by]) %>% summarise(Frequency = n()))
        if(!is.null(group.by)) {
          colnames(data2)[c(1,2)] <- c(cloneCall, group.by)
          x <- unique(data[,group.by])
          for (i in seq_along(x)) {
              sub1 <- subset(data, data[,group.by] == x[i])
              sub2 <- subset(data2, data2[,group.by] == x[i])
              merge <- merge(sub1, sub2, by=cloneCall)
              if (proportion) {
                  merge$Frequency <- merge$Frequency/length(merge$Frequency)
              }
              Con.df <- rbind.data.frame(Con.df, merge)
          }
          nsize <- Con.df %>% group_by(Con.df[,paste0(group.by, ".x")])  %>% summarise(n = n())
        } else {
          if (proportion) {
            data <- data %>%
              group_by(data[,cloneCall]) %>%
              mutate(Frequency = n()/nrow(data))
          } else {
            data <- data %>% 
              group_by(data[,cloneCall]) %>%
              mutate(Frequency = n())
          }
          Con.df <- data[,c("barcode", "CTgene", "CTnt", 
                          "CTaa", "CTstrict", "Frequency")]
          nsize <- length(Con.df)
        }
        
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
    dup <- PreMeta$barcode[which(duplicated(PreMeta$barcode))]
    PreMeta <- PreMeta[PreMeta$barcode %!in% dup,]
    barcodes <- PreMeta$barcode
    PreMeta <- PreMeta[,-1]
    rownames(PreMeta) <- barcodes
    if (group.by != "none" && addLabel) {
      location <- which(colnames(PreMeta) == "Frequency")
      colnames(PreMeta)[location] <- paste0("Frequency.", group.by)
    }
    
    warn_str <- "< 1% of barcodes match: Ensure the barcodes in 
        the Seurat object match the barcodes in the combined immune receptor
        list from scRepertoire - most common issue is the addition of the 
        prefixes corresponding to `samples` and 'ID' in the combineTCR/BCR() 
        functions"
    
    if (is_seurat_object(sc)) { 
        if (length(which(rownames(PreMeta) %in% 
                         rownames(sc[[]])))/length(rownames(sc[[]])) < 0.01) {
          warning(warn_str)
        }
        col.name <- names(PreMeta) %||% colnames(PreMeta)
        sc[[col.name]] <- PreMeta
    } else {
      rownames <- rownames(colData(sc))
      if (length(which(rownames(PreMeta) %in% 
                       rownames))/length(rownames) < 0.01) {
        warning(warn_str) }
      colData(sc) <- cbind(colData(sc), PreMeta[rownames,])[, union(colnames(colData(sc)),  colnames(PreMeta))]
      rownames(colData(sc)) <- rownames  
    }
    if (filterNA) { sc <- filteringNA(sc) }
    sc$cloneType <- factor(sc$cloneType, levels = rev(names(cloneTypes)))
    
    if(is_seurat_object(sc)) {
        sc@commands[["combineExpression"]] <- make_screp_seurat_cmd(
            call_time, sc@active.assay
        )
    }
    sc
} 











