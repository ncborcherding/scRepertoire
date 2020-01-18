#df refers to the combined contig object after the combineContig()
#seurat is the seurat object to attach
#call is the how to call the clonotype - based on genes, CDR3 nt, or CDR3 aa sequence
#groupBy is the column label in which clonotype frequency will be calculated
#cloneType are the bins for the grouping based on frequency
#' @export
combineSeurat <- function(df,
                          seurat,
                          call = c("gene", "nt", "aa"),
                          groupBy = c("none", "sample", "ID", ...),
                          cloneTypes = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded = 500)) {
    df <- if(class(df) != "list") list(df) else df
    cloneTypes <- c(None = 0, cloneTypes)
    if (length(df) == 0 | length(seurat) == 0) {
        stop("Make sure you are adding the combined contigs and Seurat object to the combineSeurat() function")
    }
    if (call == "gene") {
        call <- "CTgene"
    } else if(call == "nt") {
        call <- "CTnt"
    } else {
        call <- "CTaa"
    }
    Con.df <- NULL
    if (groupBy == "none") {
        for (i in seq_along(df)) {
            data <- data.frame(df[[i]], stringsAsFactors = F)
            data2 <- data[,c("barcode", call)]
            data2 <- unique(data2)
            data2 <- data2 %>%
                group_by(data2[,call]) %>%
                summarise(Frequency = n())
            colnames(data2)[1] <- call
            data2 <- na.omit(data2)
            data <- merge(data, data2, by = call, all = T)
            Con.df <- rbind.data.frame(Con.df, data)
            data <- NULL

        }
    }
    else if (groupBy != "none") {

        data <- bind_rows(df)
        data2 <- data[,c("barcode", call, groupBy)]
        data2 <- unique(data2)
        data2 <- data2 %>%
            group_by(data2[,call], data2[,groupBy]) %>%
            summarise(Frequency = n())
        colnames(data2)[1] <- call
        colnames(data2)[2] <- groupBy
        data2 <- na.omit(data2)
        x <- unique(data[,groupBy])
#        data2<- data2[ , -which(names(data2) %in% c(groupBy))]
        Con.df <- NULL
        for (i in seq_along(x)) {
            sub1 <- subset(data, data[,groupBy] == x[i])
            sub2 <- subset(data2, data2[,groupBy] == x[i])
            merge <- merge(sub1, sub2, by=call)
            Con.df <- rbind.data.frame(Con.df, merge)
        }
    }
    Con.df$cloneType <- NA
    for (x in seq_along(cloneTypes)) {
        names(cloneTypes)[x] <- paste0(names(cloneTypes[x]), ' (', cloneTypes[x-1], ' < X <= ', cloneTypes[x], ')')
    }
    for (i in 2:length(cloneTypes)) {
        Con.df$cloneType <- ifelse(Con.df$Frequency > cloneTypes[i-1] & Con.df$Frequency <= cloneTypes[i], names(cloneTypes[i]), Con.df$cloneType)
    }

    PreMeta <- Con.df[,c("barcode", "CTgene", "CTnt", "CTaa", "Frequency", "cloneType")]
    PreMeta <- unique(PreMeta)
    rownames(PreMeta) <- PreMeta$barcode
    seurat <- Seurat::AddMetaData(seurat, PreMeta)
}

#seurat is the seurat object to attach
#call is the how to call the clonotype - based on genes, CDR3 nt, or CDR3 aa sequence
#sequence is the specific sequence of the clonotype you would like to highlight
#' @export
highlightClonotypes <- function(seurat,
                                call = c("gene", "nt", "aa"),
                                sequence = NULL){
    if (class(seurat) != "Seurat"){
        stop("Object indicated is not of class 'Seurat', make sure you are using the correct data.")
    }
    if (call == "gene") {
        call <- "CTgene"
    } else if(call == "nt") {
        call <- "CTnt"
    } else {
        call <- "CTaa"
    }
    if (is.null(length(sequence))) {
        stop("Make sure to add  clonotype sequence(s) to me highlighted")
    }
    meta <- seurat@meta.data
    meta$highlight <- NA
    for(i in seq_along(sequence)) {
        meta$highlight <- ifelse(meta[,call] == sequence[i], paste("Clonotype", i, sep=""),  meta$highlight)
    }
    names <- rownames(meta)
    meta <- data.frame(meta[,c("highlight")])
    rownames(meta) <- names
    colnames(meta)[1] <- "highlight"
    seurat <- AddMetaData(seurat, meta)

}

#seurat is the seurat object to attach
#call is the how to call the clonotype - based on genes, CDR3 nt, or CDR3 aa sequence
#compare is the column name of the variable you'd like to compare
#facet is the additional column variable that can be used to seperate the visualization
#' @export
alluvialGraph <- function(seurat,
                          call = c("gene", "nt", "aa"),
                          compare = c("cluster", ...),
                          facet = NULL) {

    if (class(seurat) != "Seurat"){
        stop("Object indicated is not of class 'Seurat', make sure you are using the correct data.")
    }
    if (call == "gene") {
        call <- "CTgene"
    } else if(call == "nt") {
        call <- "CTnt"
    } else {
        call <- "CTaa"
    }
    if (length(compare) > 1) {
        stop("Only one comparison can be made for the sake of sanity, if you'd like to seperate by another variable, use the facet call.")
    }
    if (compare == "cluster") {
        compare <- "seurat.active.ident"
    }

    meta <- data.frame(seurat@meta.data, seurat@active.ident)

    plot <- ggplot(meta,
                   aes(axis1 = meta[,compare], axis2 = reorder(meta[,call], Frequency))) +
        geom_alluvium(aes(fill=meta[,compare])) +
        geom_stratum(fill = "grey", color = "black", lwd=0.1) +
        theme_classic() +
        labs(fill = compare) +
        ylab("Individual Clonotypes") +
        guides(fill=F) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    if (length(facet) == 1 & length(facet) < 2) {
        plot <- plot + facet_wrap(.~meta[,facet], scales="free_y")
    }
    else if (length(facet) == 0) {
        plot <- plot
    }
    else {
        stop("Only can facet on one item at a time, try multiple calls with different facets.")
    }
    return(plot)
}
