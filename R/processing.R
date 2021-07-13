#' Removing any additional prefixes to the barcodes of filtered contigs.
#'
#' @param contigs The raw loaded filtered_contig_annotation.csv
#' @param column The column in which the barcodes are listed
#' @param connector The type of character in which is attaching the defualt 
#' barcode with any other characters
#' @param num_connects The number of strings combined with the connectors
#' @examples 
#' stripBarcode(contig_list[[1]], column = 1, connector = "_", num_connects = 1)
#' @export
#' @return list with the suffixes of the barcodes removed.
stripBarcode <- function(contigs, column = 1, connector = "_", 
                            num_connects = 3) {
    count <- as.data.frame(t(data.frame(strsplit(contigs[,column], 
                            paste("['", connector, "']", sep="")), 
                            stringsAsFactors = FALSE)), 
                            stringsAsFactors = FALSE)[num_connects]
    contigs[,column] <- count
    return(contigs)
}

#' Adding variables after the combination of contigs.
#'
#' This function adds variables to the product of combineTCR() combineBCR() or 
#' expression2List() to be used in later visualizations. For each element, 
#' the function will add a column (labled by name) with the variable. 
#' The length of the variable paramater needs to match the length of 
#' the combined object.
#'
#' @examples
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' combined <- addVariable(combined, name = "batch", variables = c(1,1,1,1,2,2))
#'
#' @param df The product of combineTCR() combineBCR() or expression2List().
#' @param name The column header to add.
#' @param variables The exact values to add to each element of the list.
#' @export
#' @return list of contigs with a new column (name).
addVariable <- function(df, name = NULL, variables =  NULL) {
    if (length(df) != length(variables)) {
        stop("Make sure the variables match the length of the contig list")
    }
    for (i in seq_along(df)) {
        df[[i]][,name] <- variables[i]
    }
    return(df)
}

#' Subset the product of combineTCR() combineBCR() or expression2List()
#'
#' This function allows for the subsetting of the product of combineTCR() 
#' combineBCR() or expression2List() by the name of the individual list 
#' element. In general the names of are samples + _ + ID, allowing 
#' for users to subset the product of combineTCR(), combineBCR(), 
#' or expression2List() across a string or individual name.
#'
#' @examples
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' subset <- subsetContig(combined, name = "sample", variables = c("PX"))
#'
#' @param df The product of combineTCR(), combineBCR(), or expression2List().
#' @param name The column header you'd like to use to subset.
#' @param variables The values to subset by, must be in the names(df).

#' @export
#' @return list of contigs that have been filtered for the name parameter
subsetContig <- function(df, name, variables = NULL) {
    names2 <- NULL
    df2 <- list()
    for (i in seq_along(df)) {
        if (df[[i]][1,name] %in% variables) {
            df2 <- append(df2, list(df[[i]]))
            n2 <- names(df)[i] 
            names2 <- c(names2, n2)
        }
        else {
            next()
        }
    }
    names(df2) <- names2
    return(df2)
}

#' Allows users to take the meta data in Seurat/SCE and place it into a list 
#' that will work with all the functions
#'
#' Allows users to perform more fundamental measures of clonotype analysis 
#' using the meta data from the Seurat or SCE object. For Seurat objects the 
#' active identity is automatically added as "cluster". Remaining grouping 
#' parameters or SCE or Seurat objects must appear in the meta data.
#'
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' #Getting a sample of a Seurat object
#' screp_example <- get(data("screp_example"))
#' 
#' #Using expression2List
#' newList <- expression2List(screp_example, group = "seurat_clusters")
#' 
#' @param sc object after combineExpression().
#' @param group The column header to group the new list by
#' @importFrom stringr str_sort
#' @export
#' @return list derived from the meta data of single-cell object with 
#' elements divided by the group parameter
expression2List <- function(sc, group) {
    if (!inherits(x=sc, what ="Seurat") & 
        !inherits(x=sc, what ="SummarizedExperiment")) {
            stop("Use a Seurat or SCE object to convert into a list")
    }
    meta <- grabMeta(sc)
    unique <- str_sort(as.character(unique(meta[,group])), numeric = TRUE)
    df <- NULL
    for (i in seq_along(unique)) {
        subset <- subset(meta, meta[,group] == unique[i])
        subset <- subset(subset, !is.na(cloneType))
        df[[i]] <- subset
    }
    names(df) <- unique
    return(df)
}

#' Generate data frame to be used with circlize R package to visualize
#' clonotypes as a chord diagram. 
#' 
#' This function will take the meta data from the product of 
#' combineExpression()and generate a relational data frame to 
#' be used for a chord diagram. The output is a measure of 
#' relative clonotype overlap between groups.
#' 
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' #Getting a sample of a Seurat object
#' screp_example <- get(data("screp_example"))
#' screp_example <- combineExpression(combined, screp_example)
#' 
#' #Getting data frame output for Circilize
#' circles <- getCirclize(screp_example, groupBy = "seurat_clusters")
#' 
#' 
#' @param sc object after combineExpression().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (gene+nt).
#' @param groupBy The group header for which you would like to analyze 
#' the data.
#' @param proportion Binary will calculate relationship unique 
#' clonotypes (proportion = FALSE) or a ratio of the groupBy 
#' variable (proportion = TRUE)
#' 
#' @importFrom reshape2 dcast
#' @export
#' @return data frame of shared clonotypes between groups
#' @author Dillon Corvino, Nick Borcherding
getCirclize <- function(sc, cloneCall = "gene+nt", 
                        groupBy = NULL, proportion = FALSE) {
    meta <- grabMeta(sc)
    cloneCall <- theCall(cloneCall)
    test <- meta[, c(cloneCall, groupBy)]
    dTest <- dcast(test, test[,cloneCall] ~ test[,groupBy])
    dTest <- dTest[apply(dTest[,-1], 1, function(x) !all(x==0)),]
    dTest <- dTest[-1]
    total <- nrow(dTest)
    matrix_out <- matrix(ncol = ncol(dTest), nrow = ncol(dTest))
    for (x in seq_len(ncol(dTest))) {
        for (y in seq_len(ncol(dTest)) ){
            matrix_out[y,x] <- length(which(dTest[,x] >= 1 & dTest[,y] >= 1))
        }
    }
    colnames(matrix_out) <- colnames(dTest)
    rownames(matrix_out) <- colnames(dTest)
    
    #Need to subtract extra cells - will take the difference of the sum of the 
    #column minus and the respective cell and subtract that from the respective cell
    for (y in seq_len(ncol(matrix_out))) {
        matrix_out[y,y] <- matrix_out[y,y] - (sum(matrix_out[,y])-matrix_out[y,y])
        if (matrix_out[y,y] < 0) {
            matrix_out[y,y] <- 0
        }
    }
    output <- data.frame(from = rep(rownames(matrix_out), 
                        times = ncol(matrix_out)),
                        to = rep(colnames(matrix_out), each = nrow(matrix_out)),
                        value = as.vector(matrix_out),
                        stringsAsFactors = FALSE)
    # Reorder columns to eliminate redundant comparisons
    for (k in seq_len(nrow(output))) {
        max <- order(output[k,seq_len(2)])[1] #which is first alphabetically
        max <- output[k,max]
        min <- order(output[k,seq_len(2)])[2] #which is second alphabetically
        min <- output[k,min]
        output[k,1] <- max
        output[k,2] <- min
    }
    unique <- rownames(unique(output[,seq_len(2)])) #removing redundant comparisons
    output <- output[rownames(output) %in% unique, ]
    if (proportion == TRUE) {
        output$value <- output$value/total
    } 
    return(output)
}

