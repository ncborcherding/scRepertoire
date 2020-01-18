#df indicates the list of data frames of the contigs
#split is the point of seperation for the clonotype groupings
#call is the how to call the clonotype - based on genes, CDR3 nt, or CDR3 aa sequence
#' @export
clonalProportion <- function(df,  
                             split = c(10, 100, 1000, 10000, 30000, 100000), 
                             call = "gene") {
    df <- if(class(df) != "list") list(df) else df
    Con.df <- NULL
    if (call == "gene") {
        call <- "CTgene"
    } else if(call == "nt") {
        call <- "CTnt"
    } else {
        call <- "CTaa"
    }
    
    mat <- matrix(0, length(df), length(split), dimnames = list(names(df), paste0('[', c(1, split[-length(split)] + 1), ':', split, ']')))
    df <- lapply(df, '[[', call)
    df <- lapply(df, as.data.frame(table))
    for (i in seq_along(df)) {
        df[[i]] <- na.omit(df[[i]])
        df[[i]] <- rev(sort(as.numeric(df[[i]][,2])))
    }
    cut <- c(1, split[-length(split)] + 1)
    for (i in 1:length(split)) {
        mat[,i] <- sapply(df, function (x) sum(na.omit(x[cut[i]:split[i]])))
    }
    mat_melt <- reshape2::melt(mat)
    col <- length(unique(mat_melt$Var2))
    ggplot2::ggplot(mat_melt, aes(x=Var1, y=value, fill=Var2)) +
        geom_bar(stat = "identity", position="fill", color = "black", lwd= 0.25) +
        scale_fill_manual(name = "Clonal Indices", values = colorblind_vector(col)) +
        xlab("Samples") + 
        ylab("Occupied Repertoire Space") +
        theme_classic()
    
}