
#Ensure df is in list format
checkList <- function(df) {
    df <- if(is(df)[1] != "list") list(df) else df
    return(df)
}

#This is to check the single-cell expresison object
checkSingleObject <- function(sc) {
    if (!inherits(x=sc, what ="Seurat") | 
        inherits(x=sc, what ="SummarizedExperiment")){
        stop("Object indicated is not of class 'Seurat' or 
            'SummarizedExperiment', make sure you are using
            the correct data.") }
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

#Filtering NA contigs out of single-cell expression object
filteringNA <- function(sc) {
    meta <- grabMeta(sc)
    evalNA <- data.frame(meta[,"cloneType"])
    colnames(evalNA) <- "indicator"
    evalNA <- evalNA %>%
        transmute(indicator = ifelse(is.na(indicator), 0, 1))
    rownames(evalNA) <- rownames(meta)
    sc <- AddMetaData(sc, evalNA)
    sc <- subset(sc, cloneType != 0)
    return(sc)
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

#Calculate the Morisita Index for Overlap Analysis
morisitaIndex <- function(df, length, cloneCall, coef_matrix) {
    for (i in seq_along(length)){
        df.i <- df[[i]]
        df.i <- data.frame(table(df.i[,cloneCall]))
        colnames(df.i) <- c(cloneCall, 'Count')
        df.i[,2] <- as.numeric(df.i[,2])
        for (j in seq_along(length)){
            if (i >= j){ next }
            else { df.j <- df[[j]]
            df.j <- data.frame(table(df.j[,cloneCall]))
            colnames(df.j) <- c(cloneCall, 'Count')
            df.j[,2] <- as.numeric(df.j[,2])
            merged <- merge(df.i, df.j, by = cloneCall, all = TRUE)
            merged[is.na(merged)] <- 0
            sum.df.i <- sum(df.i[,2])
            sum.df.j <- sum(df.j[,2])
            coef.i.j <- 2 * sum(merged[,2] * merged[,3] / sum.df.j) / 
                sum.df.j/((sum((df.i[,2] / sum.df.i)^2) + 
                sum((df.j[,2] / sum.df.j)^2)))
            coef_matrix[i,j] <- coef.i.j } } }
    return(coef_matrix)
}

#Calculate the Overlap Coefficient for Overlap Analysis
overlapIndex <- function(df, length, cloneCall, coef_matrix) {
    for (i in seq_along(length)){
        df.i <- df[[i]]
        df.i <- df.i[,c("barcode",cloneCall)]
        df.i_unique <- df.i[!duplicated(df.i$barcode),]
        for (j in seq_along(length)){
            if (i >= j){ next }
            else { df.j <- df[[j]]
            df.j <- df.j[,c("barcode",cloneCall)]
            df.j_unique <- df.j[!duplicated(df.j$barcode),]
            overlap <- length(intersect(df.i_unique[,cloneCall], 
                                        df.j_unique[,cloneCall]))
            coef_matrix[i,j] <- 
                overlap/min(length(df.i_unique[,cloneCall]), 
                length(df.j_unique[,cloneCall])) } } }
    return(coef_matrix)
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