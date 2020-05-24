#' Removing any additional prefixes to the barcodes of filtered contigs.
#'
#' @param contigs The raw loaded filtered_contig_annotation.csv
#' @param column The column in which the barcodes are listed
#' @param connector The type of character in which is attaching the defualt 
#' barcode with any other characters
#' @param num_connects The number of strings combined with the connectors
#' @examples 
#' stripBarcode(contig_list[[1]], column = 1, connector = "_", num_connects = 3)
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
#' This function adds variables to the product of combineContig() to be used in 
#' later visualizations. For each element in the combinContig(), the function 
#' will add a column (labled by name) with the variable. The length of the
#' variable paramater needs to match the length of the combineContig() object.
#'
#' @examples
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' combined <- addVariable(combined, name = "batch", variables = c(1,1,1,1,2,2))
#'
#' @param df The product of CombineContig().
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

#' Subset the combineContig() product.
#'
#' This function allows for the subsetting of the product of combineContig() by 
#' the name of the individual list element. In general the names of are samples 
#' + _ + ID, allowing for users to subset the product of combineContig() across 
#' a string or individual name.
#'
#' @examples
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' subset <- subsetContig(combined, name = "sample", variables = c("PX"))
#'
#' @param df The product of CombineContig().
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

# This suppressing outputs for using dput()
quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}

#This is to help sort the type of clonotype data to use
theCall <- function(x) {
if (x == "gene") {
    x <- "CTgene"
} else if(x == "nt") {
    x <- "CTnt"
} else if (x == "aa") {
    x <- "CTaa"
} else if (x == "gene+nt") {
    x <- "CTstrict"
}
    return(x)
}

#' Allows users to take the meta data in seurat/SCE and place it into a list 
#' that will work with all the functions
#'
#' Allows users to perform more fundamental measures of clonotype analysis 
#' using the meta data from the seurat or SCE object. For Seurat objects the 
#' active identity is automatically added as "cluster". Reamining grouping 
#' parameters or SCE or Seurat objects must appear in the meta data.
#'
#' @examples
#' \donttest{
#' newList <- expression2List(seurat, group = "cluster")
#' }
#' @param sc object after combineExpression().
#' @param group The column header to group the new list by
#' @importFrom stringr str_sort
#' @export
#' @return list derived from the meta data of single-cell object with 
#' elements divided by the group parameter
expression2List <- function(sc, group) {
    if (!inherits(x=sc, what ="Seurat") & 
        !inherits(x=sc, what ="SummarizedExperiment")) {
            stop("Use a seurat or SCE object to convert into a list")
    }
    meta <- grabMeta(sc)
    unique <- str_sort(as.character(unique(meta[,group])), numeric = TRUE)
    df <- NULL
    for (i in seq_along(unique)) {
        subset <- subset(meta, meta[,"cluster"] == unique[i])
        df[[i]] <- subset
    }
    names(df) <- unique
    return(df)
}

#This is to check the single-cell expresison object
checkSingleObject <- function(sc) {
    if (!inherits(x=sc, what ="Seurat") | 
        inherits(x=sc, what ="SummarizedExperiment")){
    stop("Object indicated is not of class 'Seurat' or 
        'SummarizedExperiment', make sure you are using the correct data.") }
}

#This is to grab the meta data from a seurat or SCE object
grabMeta <- function(sc) {
    if (inherits(x=sc, what ="Seurat")) {
        meta <- data.frame(sc[[]], Idents(sc))
        colnames(meta)[length(meta)] <- "cluster"
    }
    else if (inherits(x=sc, what ="SummarizedExperiment")){
        meta <- sc@metadata[[1]]
    }
    return(meta)
}

#This is to add the sample and ID prefixes for combineBCR()/TCR()
modifyBarcodes <- function(df, samples, ID) {
    out <- NULL
    for (x in seq_along(df)) {
        data <- df[[x]]
        data$barcode <- paste0(samples[x], "_", ID[x], "_", data$barcode)
        out[[x]] <- data }
    return(out)
}

#Removing barcodes with NA recovered
removingNA <- function(final) {
    for(i in seq_along(final)) {
        final[[i]] <- na.omit(final[[i]])}
    return(final)
}

#Removing barcodes with > 2 clones recovered
removingMulti <- function(final){
    for(i in seq_along(final)) {
        final[[i]] <- filter(final[[i]], !grepl(";",CTnt))}
    return(final)
}

#Check the format of the cell barcode inputs and parameter lengthsd
checkContigBarcodes <- function(df, samples, ID) {
    count <- length(unlist(strsplit(df[[1]]$barcode[1], "[-]")))
    count2 <- length(unlist(strsplit(df[[1]]$barcode[1], "[_]")))
    if (count > 2 | count2 > 2) {
        stop("Seems to be an error in the naming of the contigs, ensure 
            the barcodes are labeled like, AAACGGGAGATGGCGT-1 or 
            AAACGGGAGATGGCGT, use stripBarcode to get the basic 
            format", call. = FALSE)
    } else if (length(df) != length(samples) | length(df) != length(ID)) {
        stop("Make sure the sample and ID labels match the length of the 
            list of data frames (df).", call. = FALSE) }
}

#Caclulating diversity using Vegan R package
#' @importFrom vegan diversity estimateR
diversityCall <- function(data) {
    w <- diversity(data[,"Freq"], index = "shannon")
    x <- diversity(data[,"Freq"], index = "invsimpson")
    y <- estimateR(data[,"Freq"])[2] #Chao
    z <- estimateR(data[,"Freq"])[4] #ACE
    out <- c(w,x,y,z)
    return(out)
}

#Organizing list of contigs for vizualization
parseContigs <- function(df, i, names, cloneCall) {
    data <- df[[i]]
    data1 <- data %>% group_by(data[,cloneCall]) %>%
        summarise(Abundance=n())
    colnames(data1)[1] <- cloneCall
    data1$values <- names[i]
    return(data1)
}
