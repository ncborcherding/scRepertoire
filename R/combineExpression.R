#' Adding clone information to a single-cell object
#'
#' This function adds the immune receptor information to the Seurat or 
#' SCE object to the meta data. By default this function also calculates 
#' the frequencies and proportion of the clones by sequencing 
#' run (\strong{group.by} = NULL). To change how the frequencies/proportions
#' are calculated, select a column header for the \strong{group.by} variable. 
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
#' 
#' #Using combineExpresion()
#' scRep_example <- combineExpression(combined, scRep_example)
#' 
#' @param input.data The product of \code{\link{combineTCR}}, 
#' \code{\link{combineBCR}} or a list of 
#' both c(\code{\link{combineTCR}}, \code{\link{combineBCR}}).
#' @param sc.data The Seurat or Single-Cell Experiment (SCE) object to attach
#' @param cloneCall How to call the clone - VDJC gene (\strong{gene}), 
#' CDR3 nucleotide (\strong{nt}), CDR3 amino acid (\strong{aa}),
#' VDJC gene + CDR3 nucleotide (\strong{strict}) or a custom variable 
#' in the data. 
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL".
#' @param group.by The column label in the combined clones in which 
#' clone frequency will be calculated. \strong{NULL} or \strong{"none"} will 
#' keep the format of input.data.
#' @param proportion Whether to proportion (\strong{TRUE}) or total 
#' frequency (\strong{FALSE}) of the clone based on the group.by variable. 
#' @param cloneSize The bins for the grouping based on proportion or frequency. 
#' If proportion is \strong{FALSE} and the cloneSizes are not set high enough
#' based on frequency, the upper limit of cloneSizes will be automatically
#' updated.S
#' @param filterNA Method to subset Seurat/SCE object of barcodes without 
#' clone information
#' @param addLabel This will add a label to the frequency header, allowing
#' the user to try multiple group.by variables or recalculate frequencies after 
#' subsetting the data.
#' @importFrom dplyr bind_rows %>% summarise left_join mutate select n all_of coalesce
#' @importFrom  rlang %||% sym :=
#' @importFrom SummarizedExperiment colData<- colData
#' @importFrom S4Vectors DataFrame
#' @export
#' @concept SC_Functions
#' @return Single-cell object with clone information added to meta data
#' information
#' 
combineExpression <- function(input.data, 
                              sc.data, 
                              cloneCall ="strict", 
                              chain = "both", 
                              group.by = NULL, 
                              proportion = TRUE, 
                              filterNA = FALSE,
                              cloneSize = c(Rare = 1e-4,Small = 0.001,Medium = 0.01,Large = 0.1,Hyperexpanded = 1),
                              addLabel = FALSE) {
    call_time <- Sys.time()
  
    options( dplyr.summarise.inform = FALSE )
    if (!proportion && any(cloneSize < 1)) {
        stop("Adjust the cloneSize parameter - there are groupings < 1")
    }
    cloneSize <- c(None = 0, cloneSize)
    
    cloneCall <- .theCall(input.data, cloneCall)
    if (chain != "both") {
      #Retain the full clone information
      full.clone <- lapply(input.data, function(x) {
                        x[,c("barcode", cloneCall)]
      full.clone <- bind_rows(full.clone)
      })
      for(i in seq_along(input.data)) {
        input.data[[i]] <- .off.the.chain(input.data[[i]], chain, cloneCall)
      }
    }
    input.data <- .checkList(input.data)
    
    #Getting Summaries of clones from combineTCR() or combineBCR()
    Con.df <- NULL
    meta <- .grabMeta(sc.data)
    cell.names <- rownames(meta)
    if (is.null(group.by) || group.by == "none") {
        for (i in seq_along(input.data)) {
      
            data <- data.frame(input.data[[i]], stringsAsFactors = FALSE)
            data2 <- unique(data[,c("barcode", cloneCall)])
            #This ensures all calculations are based on the cells in the SCO
            data2 <- na.omit(data2[data2[,"barcode"] %in% cell.names,])
            data2 <- data2 %>% 
                        group_by(data2[,cloneCall]) %>%
                        summarise(clonalProportion = dplyr::n()/nrow(data2), 
                                  clonalFrequency = dplyr::n())
            colnames(data2)[1] <- cloneCall
            data <- merge(data, data2, by = cloneCall, all = TRUE)
            if ( cloneCall %!in% c("CTgene", "CTnt", "CTaa", "CTstrict") ) {
              data <- data[,c("barcode", "CTgene", "CTnt",
                              "CTaa", "CTstrict", cloneCall,
                              "clonalProportion", "clonalFrequency")]
            } else {
              data <- data[,c("barcode", "CTgene", "CTnt", 
                              "CTaa", "CTstrict",
                              "clonalProportion", "clonalFrequency")] }
            Con.df <- rbind.data.frame(Con.df, data)
        }
    } else if (group.by != "none" || !is.null(group.by)) {
        data <- data.frame(bind_rows(input.data), stringsAsFactors = FALSE)
        data2 <- na.omit(unique(data[,c("barcode", cloneCall, group.by)]))
        #This ensures all calculations are based on the cells in the SCO
        data2 <- data2[data2[,"barcode"] %in% cell.names, ]
        data2 <- as.data.frame(data2 %>% 
                                  group_by(data2[,cloneCall], data2[,group.by]) %>% 
                                  summarise(clonalProportion = dplyr::n()/nrow(data2), 
                                            clonalFrequency = dplyr::n())
        )
        
        colnames(data2)[c(1,2)] <- c(cloneCall, group.by)
        data <- merge(data, data2, by = c(cloneCall, group.by), all = TRUE)
        if ( cloneCall %!in% c("CTgene", "CTnt", "CTaa", "CTstrict") ) {
              Con.df <- data[,c("barcode", "CTgene", "CTnt",
                              "CTaa", "CTstrict", cloneCall,
                              "clonalProportion", "clonalFrequency")]
            } else {
              Con.df <- data[,c("barcode", "CTgene", "CTnt", 
                              "CTaa", "CTstrict",
                              "clonalProportion", "clonalFrequency")] }
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
    if ( cloneCall %!in% c("CTgene", "CTnt", 
                         "CTaa", "CTstrict") ) {
      PreMeta <- unique(Con.df[,c("barcode", "CTgene", "CTnt", 
                                  "CTaa", "CTstrict", cloneCall, 
                                  "clonalProportion", "clonalFrequency", "cloneSize")])
    } else {
      PreMeta <- unique(Con.df[,c("barcode", "CTgene", "CTnt", 
                                "CTaa", "CTstrict", "clonalProportion", 
                                "clonalFrequency", "cloneSize")])
    }
    #Removing any duplicate barcodes, should not be an issue
    dup <- PreMeta$barcode[which(duplicated(PreMeta$barcode))]
    PreMeta <- PreMeta[PreMeta$barcode %!in% dup,]
    
    #Re-adding full clones 
    if (chain != "both") {
      clone_sym <- sym(cloneCall)
      PreMeta <- PreMeta %>%
        left_join(full.clone, by = "barcode", suffix = c("", ".from_full_clones")) %>%
        mutate(!!column_sym := coalesce(!!sym(paste0(cloneCall, ".from_full_clones")), !!column_sym)) %>%
        select(-all_of(paste0(cloneCall, ".from_full_clones")))
    }
    barcodes <- PreMeta$barcode
    PreMeta <- PreMeta[,-1]
    rownames(PreMeta) <- barcodes
    if (group.by != "none" && addLabel) {
      location <- which(colnames(PreMeta) %in% c("clonalProportion", 
                          "clonalFrequency"))
      colnames(PreMeta)[location] <- paste0(c("clonalProportion", 
                                            "clonalFrequency"), group.by)
    }
    
    if (is_seurat_object(sc.data)) { 
        if (length(which(rownames(PreMeta) %in% 
                         rownames(sc.data[[]])))/length(rownames(sc.data[[]])) < 0.01) {
          warning(.warn_str)
        }
        col.name <- names(PreMeta) %||% colnames(PreMeta)
        sc.data[[col.name]] <- PreMeta
    } else {
      rownames <- rownames(colData(sc.data))
      if (length(which(rownames(PreMeta) %in% 
                       rownames))/length(rownames) < 0.01) {
        warning(.warn_str) }
      
      combined_col_names <- unique(c(colnames(colData(sc.data)), colnames(PreMeta)))
      full_data <- merge(colData(sc.data), PreMeta[rownames, , drop = FALSE], by = "row.names", all.x = TRUE)
      rownames(full_data) <- full_data[, 1]
      full_data  <- full_data[, -1]
      colData(sc.data) <- DataFrame(full_data[, combined_col_names])
      
      rownames(colData(sc.data)) <- rownames  
    }
    if (filterNA) { 
      sc.data <- .filteringNA(sc.data) 
    }
    sc.data$cloneSize <- factor(sc.data$cloneSize, levels = rev(names(cloneSize)))
    
    if(is_seurat_object(sc.data)) {
        sc.data@commands[["combineExpression"]] <- make_screp_seurat_cmd(
            call_time, sc.data@active.assay
        )
    }
    return(sc.data)
} 


.warn_str <- "< 1% of barcodes match: Ensure the barcodes in the single-cell object match the barcodes in the combined immune receptor output from scRepertoire. If getting this error, please check https://www.borch.dev/uploads/screpertoire/articles/faq."









