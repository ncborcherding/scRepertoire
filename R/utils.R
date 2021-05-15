
#Remove list elements that contain all NA values
checkBlanks <- function(df, cloneCall) {
    for (i in seq_along(df)) {
        if (length(df[[i]][,cloneCall]) == length(which(is.na(df[[i]][,cloneCall])))) {
            df[[i]] <- NULL
        } else {
            next()
        }
    }
    return(df)
}


#Ensure df is in list format
checkList <- function(df) {
    df <- if(is(df)[1] != "list") list(df) else df
    return(df)
}

#This is to check the single-cell expresison object
checkSingleObject <- function(sc) {
    if (!inherits(x=sc, what ="Seurat") & 
        !inherits(x=sc, what ="SummarizedExperiment")){
        stop("Object indicated is not of class 'Seurat' or 
            'SummarizedExperiment', make sure you are using
            the correct data.") }
    }

#This is to grab the meta data from a seurat or SCE object
#' @importFrom SummarizedExperiment colData 
grabMeta <- function(sc) {
    if (inherits(x=sc, what ="Seurat")) {
        meta <- data.frame(sc[[]], slot(sc, "active.ident"))
        colnames(meta)[length(meta)] <- "cluster"
    }
    else if (inherits(x=sc, what ="SummarizedExperiment")){
        meta <- data.frame(colData(sc))
        rownames(meta) <- sc@colData@rownames
        clu <- which(colnames(meta) == "ident")
        colnames(meta)[clu] <- "cluster"
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

#Removing extra clonotypes in barcodes with > 2
#' @import dplyr
filteringMulti <- function(x) {
    table <- subset(as.data.frame(table(x$barcode)), Freq > 2)
    barcodes <- as.character(unique(table$Var1))
    multichain <- NULL
    for (j in seq_along(barcodes)) {
        chain <- x[x$barcode == barcodes[j],] %>% 
            group_by(barcode, chain) %>% top_n(n = 1, wt = reads)
        multichain <- rbind(multichain, chain) }
    `%!in%` = Negate(`%in%`)
    x <- subset(x, barcode %!in% barcodes)
    x <- rbind(x, multichain) 
    return(x)
}



#Filtering NA contigs out of single-cell expression object
#' @import dplyr
filteringNA <- function(sc) {
    meta <- grabMeta(sc)
    evalNA <- data.frame(meta[,"cloneType"])
    colnames(evalNA) <- "indicator"
    evalNA <- evalNA %>%
        transmute(indicator = ifelse(is.na(indicator), 0, 1))
    rownames(evalNA) <- rownames(meta)
    if (inherits(x=sc, what ="cell_data_set")){
      colData(sc)[["evalNA"]]<-evalNA
      return(sc[, !is.na(sc$cloneType)])
    }else{
      col.name <- names(evalNA) %||% colnames(evalNA)
      sc[[col.name]] <- evalNA
      sc <- subset(sc, cloneType != 0)
      return(sc)
    }
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
#' @author Massimo Andreatta, Nick Borcherding
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
            X <- sum(merged[,2])
            Y <- sum(merged[,3])
            sum.df.i <- sum(df.i[,2]^2)
            sum.df.j <- sum(df.j[,2]^2)
            
            num <- 2 * sum(merged[, 2] * merged[, 3])
            den <- ((sum.df.i / (X^2) + sum.df.j / (Y^2)) * X * Y)
                
            coef.i.j <- num/den
            coef_matrix[i,j] <- coef.i.j
            }
        }
    }
    return(coef_matrix)
}

#Calculate the Overlap Coefficient for Overlap Analysis
#' @author Nick Bormann, Nick Borcherding
overlapIndex <- function(df, length, cloneCall, coef_matrix) {
    for (i in seq_along(length)){
        df.i <- df[[i]]
        df.i <- df.i[,c("barcode",cloneCall)]
        df.i_unique <- df.i[!duplicated(df.i[,cloneCall]),]
        for (j in seq_along(length)){
            if (i >= j){ next }
            else { df.j <- df[[j]]
            df.j <- df.j[,c("barcode",cloneCall)]
            df.j_unique <- df.j[!duplicated(df.j[,cloneCall]),]
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

# This is to help sort the type of clonotype data to use
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

# Assigning positions for TCR contig data
#' @author Gloria Kraus, Nick Bormann, Nick Borcherding
parseTCR <- function(Con.df, unique_df, data2) {
    for (y in seq_along(unique_df)){
        barcode.i <- Con.df$barcode[y]
        location.i <- which(barcode.i == data2$barcode)
        if (length(location.i) == 2){
            if (is.na(data2[location.i[1],c("TCR1")])) {
                Con.df[y,tcr2_lines]<-data2[location.i[1],data2_lines]
                if (is.na(data2[location.i[2],c("TCR2")])) {
                  Con.df[y,tcr1_lines]<-data2[location.i[2],data1_lines]
                }
            } else {
              if(!is.na(data2[location.i[1],c("TCR1")])) {
                Con.df[y,tcr1_lines]<-data2[location.i[1],data1_lines]
              } 
              if (!is.na(data2[location.i[2],c("TCR2")])) {
                Con.df[y,tcr2_lines]<-data2[location.i[2],data2_lines] 
              }
            }
        } else if (length(location.i) == 3) { 
            if (is.na(data2[location.i[1],c("TCR1")])) { 
                Con.df[y,tcr2_lines]<-data2[location.i[1],data2_lines] 
                if (is.na(data2[location.i[2],c("TCR1")])) { 
                    TRdf <- paste(Con.df[y, tcr2_lines],
                        data2[location.i[2], data2_lines],sep=";") 
                    Con.df[y,tcr2_lines] <- TRdf 
                    Con.df[y,tcr1_lines] <- data2[location.i[3],data1_lines] 
                } else { # if the 2nd location is occupied by TRA
                    Con.df[y,tcr1_lines] <- data2[location.i[2],data1_lines] 
                    if (is.na(data2[location.i[3],c("TCR1")])) { 
                        TRdf <- paste(Con.df[y, tcr2_lines],
                            data2[location.i[3], data2_lines],sep=";") 
                        Con.df[y,tcr2_lines] <- TRdf 
                    } else { # if the 3rd location is occupied by TRA
                        TRdf <- paste(Con.df[y, tcr1_lines],
                            data2[location.i[3], data1_lines],sep=";") 
                        Con.df[y,tcr1_lines] <- TRdf }}
            } else { # if 1st location is occupied by TRA
                Con.df[y,tcr1_lines] <- data2[location.i[1],data1_lines] 
                if (is.na(data2[location.i[2],c("TCR1")])) { 
                    if (is.na(data2[location.i[3],c("TCR1")])) { #Two TRB chains
                        TRdf <- paste(data2[location.i[2], data2_lines],
                            data2[location.i[3], data2_lines],sep=";") 
                        Con.df[y,tcr2_lines] <- TRdf 
                    } else if (!is.na(data2[location.i[3],c("TCR1")])) { # if TRA is on 3rd location
                        TRdf <- paste(Con.df[y, tcr1_lines],
                            data2[location.i[3],data1_lines],sep=";") 
                        Con.df[y,tcr1_lines] <- TRdf 
                        Con.df[y,tcr2_lines] <- data2[location.i[2],data2_lines] 
                        }
                } else { # if TRA is on 2nd location
                    TRdf <- paste(Con.df[y, tcr1_lines],
                        data2[location.i[2], data1_lines],sep=";") 
                    Con.df[y,tcr1_lines] <- TRdf 
                    Con.df[y,tcr2_lines] <- data2[location.i[3],data2_lines]}}
        } else if (length(location.i) == 1) {
            chain.i <- data2$chain[location.i]
            if (chain.i == "TRA"){
                Con.df[y,tcr1_lines] <- data2[location.i[1],data1_lines]
            } else {Con.df[y,tcr2_lines] <- data2[location.i[1],data2_lines]}}}
return(Con.df)}

#Assiging positions for BCR contig data
#' @author Gloria Kraus, Nick Bormann, Nick Borcherding
parseBCR <- function(Con.df, unique_df, data2) {
    for (y in seq_along(unique_df)){
        barcode.i <- Con.df$barcode[y]
        location.i <- which(barcode.i == data2$barcode)
        
        if (length(location.i) == 2){
            if (is.na(data2[location.i[1],c("IGHct")])) {
                Con.df[y,heavy_lines]<-data2[location.i[2], h_lines]
                if (is.na(data2[location.i[1],c("IGKct")])) {
                    Con.df[y,light_lines]<-data2[location.i[1], l_lines]
                } else {
                    Con.df[y,light_lines]<-data2[location.i[1], k_lines]}
            } else { 
              Con.df[y,heavy_lines]<-data2[location.i[1], h_lines]
              if (!is.na(data2[location.i[2],c("IGKct")])) {
                Con.df[y,light_lines]<-data2[location.i[2],k_lines]
              } else {
                Con.df[y,light_lines]<- data2[location.i[2],l_lines]
                }}
          
        } else if (length(location.i) == 1) {
            chain.i <- data2$chain[location.i]
            if (chain.i == "IGH"){
                Con.df[y,heavy_lines]<-data2[location.i[1],h_lines]
            } else if (chain.i == "IGL") {
                Con.df[y,light_lines]<- data2[location.i[1],l_lines]}
            else {
                Con.df[y,light_lines]<-data2[location.i[1], k_lines]
            }}}
    return(Con.df)
}

#Assign T/B cell chains and celltypes for combineTCR() and lengthContig
cellT <- function(cells) {
    if (cells == "T-AB") { 
        chain1 <- "TRA"
        chain2 <- "TRB" 
        cellType <- "T-AB" 
    } else if (cells == "T-GD") {
        chain1 <- "TRD"
        chain2 <- "TRG"
        cellType <- "T-GD" 
    } else if (cells == "B") {
        chain1 <- "IGH"
        chain2 <- "IGL"
        cellType <- "B" 
    }
    return(list(chain1, chain2, cellType))
}


#Producing a data frame to visualize for lengthContig()
lengthDF <- function(df, cloneCall, chains, group, c1, c2){
    Con.df <- NULL
    names <- names(df)
    if (chains == "combined") {
            for (i in seq_along(df)) {
                length <- nchar(gsub("_", "", df[[i]][,cloneCall]))
                val <- df[[i]][,cloneCall]
                if (!is.null(group)) { 
                    cols <- df[[i]][,group]
                    data <- na.omit(data.frame(length, val, cols, names[i]))
                    colnames(data) <- c("length", "CT", group, "values")
                    Con.df<- rbind.data.frame(Con.df, data) 
                } else {
                    data <- na.omit(data.frame(length, val, names[i]))
                    colnames(data) <- c("length", "CT", "values")
                    Con.df<- rbind.data.frame(Con.df, data) }}
    } else if (chains == "single") {
            for (x in seq_along(df)) {
                strings <- df[[x]][,cloneCall]
                strings <- as.data.frame(str_split(strings, "_", 
                        simplify = TRUE), stringsAsFactors = FALSE)
                val1 <- strings[,1]
                for (i in seq_along(val1)) {
                    if (grepl(";", val1[i]) == TRUE) {
                        val1[i] <- str_split(val1, ";", simplify = TRUE)[1] 
                    } else { next() } }
                val2 <- strings[,2]
                for (i in seq_along(val2)) {
                    if (grepl(";", val2[i]) == TRUE) {
                        val2[i] <- str_split(val2, ";", simplify = TRUE)[1]
                    } else { next() } }
                chain1 <- nchar(val1)
                chain2 <- nchar(val2)
                if (!is.null(group)) {
                    cols1 <- df[[x]][,group]
                    data1 <- data.frame(chain1, val1, names[x], c1, cols1)
                    colnames(data1)<-c("length","CT","values","chain",group)
                    cols2 <- df[[x]][,group]
                    data2 <- data.frame(chain2, val2, names[x], c2, cols2)
                    colnames(data2)<-c("length","CT","values","chain",group)
                }else if (is.null(group)){
                    data1 <- data.frame(chain1, val1, names[x], c1)
                    colnames(data1) <- c("length", "CT", "values", "chain")
                    data2 <- data.frame(chain2, val2, names[x], c2)
                    colnames(data2) <- c("length", "CT", "values", "chain")}
                data <- na.omit(rbind(data1, data2))
                data <- subset(data, CT != "NA" | CT != "")
                Con.df<- rbind.data.frame(Con.df, data) }}
return(Con.df)}

#General combination of nucleotide, aa, and gene sequences for T/B cells
assignCT <- function(cellType, Con.df) {
    if (cellType %in% c("T-AB", "T-GD")) {
        Con.df$CTgene <- paste(Con.df$TCR1, Con.df$TCR2, sep="_")
        Con.df$CTnt <- paste(Con.df$cdr3_nt1, Con.df$cdr3_nt2, sep="_")
        Con.df$CTaa <- paste(Con.df$cdr3_aa1, Con.df$cdr3_aa2, sep="_")
        Con.df$CTstrict <- paste(Con.df$TCR1, Con.df$cdr3_nt1, 
            Con.df$TCR2, Con.df$cdr3_nt2, sep="_")
    } else {
        Con.df$CTgene <- paste(Con.df$IGH, Con.df$IGLC, sep="_")
        Con.df$CTnt <- paste(Con.df$cdr3_nt1, Con.df$cdr3_nt2, sep="_")
        Con.df$CTaa <- paste(Con.df$cdr3_aa1, Con.df$cdr3_aa2, sep="_") }
return(Con.df)
}

#Sorting the V/D/J/C gene sequences for T and B cells
makeGenes <- function(cellType, data2, chain1, chain2) {
    if(cellType %in% c("T-AB", "T-GD")) {
        data2 <- data2 %>% 
            mutate(TCR1 = ifelse(chain == chain1, paste(with(data2, 
            interaction(v_gene,  j_gene, c_gene))), NA)) %>%
            mutate(TCR2 = ifelse(chain == chain2, paste(with(data2, 
            interaction(v_gene,  j_gene, d_gene, c_gene))), NA))
    }
    else {
        data2 <- data2 %>% 
            mutate(IGKct = ifelse(chain == "IGK", paste(with(data2, 
            interaction(v_gene,  j_gene, c_gene))), NA)) %>%
            mutate(IGLct = ifelse(chain == "IGL", paste(with(data2, 
            interaction(v_gene,  j_gene, c_gene))), NA)) %>%
            mutate(IGHct = ifelse(chain == "IGH", paste(with(data2, 
            interaction(v_gene, j_gene, d_gene, c_gene))), NA))
    }
    return(data2)
    
}
