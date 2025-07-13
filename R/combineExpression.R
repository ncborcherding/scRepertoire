#' Adding Clonal Information to Single-Cell Object
#'
#' This function adds the immune receptor information to the Seurat or 
#' SCE object to the meta data. By default this function also calculates 
#' the frequencies and proportion of the clones by sequencing 
#' run (`group.by` = NULL). To change how the frequencies/proportions
#' are calculated, select a column header for the `group.by` variable. 
#' Importantly, before using [combineExpression()] ensure the 
#' barcodes of the single-cell object object match the barcodes in the output 
#' of the [combineTCR()] or [combineBCR()]. 
#'
#' @examples
#' # Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' # Getting a sample of a Seurat object
#' scRep_example <- get(data("scRep_example"))
#' 
#' # Using combineExpresion()
#' scRep_example <- combineExpression(combined, scRep_example)
#' 
#' @param input.data The product of [combineTCR()], 
#' [combineBCR()] or a list of 
#' both c([combineTCR()], [combineBCR()]).
#' @param sc.data The Seurat or Single-Cell Experiment (SCE) object to attach
#' @param cloneCall Defines the clonal sequence grouping. Accepted values 
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino 
#' acid sequence), or `strict` (VDJC). A custom column header can also be used.
#' @param chain The TCR/BCR chain to use. Use `both` to include both chains 
#' (e.g., TRA/TRB). Accepted values: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGL` 
#' (for both light chains), `both`.
#' @param group.by A column header in lists to group the analysis 
#' by (e.g., "sample", "treatment"). If `NULL`, will be based on the list element.
#' @param proportion Whether to proportion (**TRUE**) or total 
#' frequency (**FALSE**) of the clone based on the group.by variable. 
#' @param cloneSize The bins for the grouping based on proportion or frequency. 
#' If proportion is **FALSE** and the cloneSizes are not set high enough
#' based on frequency, the upper limit of cloneSizes will be automatically
#' updated.S
#' @param filterNA Method to subset Seurat/SCE object of barcodes without 
#' clone information
#' @param addLabel This will add a label to the frequency header, allowing
#' the user to try multiple group.by variables or recalculate frequencies after 
#' subsetting the data.
#' @importFrom dplyr left_join all_of coalesce
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
        x[, c("barcode", cloneCall)]
      })
      full.clone <- bind_rows(full.clone)
      for(i in seq_along(input.data)) {
        input.data[[i]] <- .offTheChain(input.data[[i]], chain, cloneCall)
      }
    }
    input.data <- .checkList(input.data)
    
    #Getting Summaries of clones from combineTCR() or combineBCR()
    Con.df <- NULL
    meta <- .grabMeta(sc.data)
    cell.names <- rownames(meta)

    conDfColnamesNoCloneSize <- unique(c(
        "barcode", CT_lines, cloneCall, "clonalProportion", "clonalFrequency"
    ))

    # Computes the clonalProportion and clonalFrequency for each clone
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
            data <- data[, conDfColnamesNoCloneSize]
            Con.df <- rbind.data.frame(Con.df, data)
        }

    } else {
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
        Con.df <- data[, conDfColnamesNoCloneSize]
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

    cloneRatioColname <- ifelse(proportion, "clonalProportion", "clonalFrequency")

    #Assigning cloneSize
    for (i in 2:length(cloneSize)) { 
        Con.df$cloneSize <- ifelse(Con.df[, cloneRatioColname] > cloneSize[i-1] & 
                                   Con.df[, cloneRatioColname] <= cloneSize[i], 
                                   names(cloneSize[i]), 
                                   Con.df$cloneSize)
    }
    
    #Formating the meta data to add and removing any duplicate barcodes
    PreMeta <- unique(Con.df[, c(conDfColnamesNoCloneSize, "cloneSize")])
    dup <- PreMeta$barcode[which(duplicated(PreMeta$barcode))]
    PreMeta <- PreMeta[!PreMeta$barcode %in% dup,]
    
    #Re-adding full clones 
    if (chain != "both") {
      clone_sym <- sym(cloneCall)
      PreMeta <- PreMeta %>%
        left_join(full.clone, by = "barcode", suffix = c("", ".from_full_clones")) %>%
        mutate(!!clone_sym := coalesce(!!sym(paste0(cloneCall, ".from_full_clones")), !!clone_sym)) %>%
        dplyr::select(-all_of(paste0(cloneCall, ".from_full_clones")))
    }
    barcodes <- PreMeta$barcode
    PreMeta <- PreMeta[,-1]
    rownames(PreMeta) <- barcodes
    if (!is.null(group.by) && group.by != "none" && addLabel) {
      location <- which(colnames(PreMeta) %in% c("clonalProportion", 
                          "clonalFrequency"))
      colnames(PreMeta)[location] <- paste0(c("clonalProportion", 
                                            "clonalFrequency"), group.by)
    }
    
    if (.is.seurat.object(sc.data)) { 
        if (length(which(rownames(PreMeta) %in% 
                         rownames(sc.data[[]])))/length(rownames(sc.data[[]])) < 0.01) {
          warning(getHighBarcodeMismatchWarning())
        }
        col.name <- names(PreMeta) %||% colnames(PreMeta)
        sc.data[[col.name]] <- PreMeta
    } else {
      rownames <- rownames(colData(sc.data))
      if (length(which(rownames(PreMeta) %in% 
                       rownames))/length(rownames) < 0.01) {
        warning(getHighBarcodeMismatchWarning()) }
      
      combined_col_names <- unique(c(colnames(colData(sc.data)), colnames(PreMeta)))
      full_data <- merge(colData(sc.data), PreMeta[rownames, , drop = FALSE], by = "row.names", all.x = TRUE)
      # at this point, the rows in full_data are shuffled. match back with the original colData
      full_data <- full_data[match(rownames, full_data[,1]), ]
      rownames(full_data) <- full_data[, 1]
      full_data  <- full_data[, -1]
      colData(sc.data) <- DataFrame(full_data[, combined_col_names])  
    }
    if (filterNA) { 
      sc.data <- .filteringNA(sc.data) 
    }
    sc.data$cloneSize <- factor(sc.data$cloneSize, levels = rev(names(cloneSize)))
    
    if(.is.seurat.object(sc.data)) {
        sc.data@commands[["combineExpression"]] <- .makeScrepSeurat(
              call_time, sc.data@active.assay)
    }
    return(sc.data)
}

getHighBarcodeMismatchWarning <- function() paste(
    "< 1% of barcodes match: Ensure the barcodes in the single-cell object",
    "match the barcodes in the combined immune receptor output from",
    "scRepertoire. If getting this error, please check",
    "https://www.borch.dev/uploads/screpertoire/articles/faq"
)
