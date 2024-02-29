# readability functions
"%!in%" <- Negate("%in%")
is_seurat_object <- function(obj) inherits(obj, "Seurat")
is_se_object <- function(obj) inherits(obj, "SummarizedExperiment")
is_seurat_or_se_object <- function(obj) {
    is_seurat_object(obj) || is_se_object(obj)
}

#Use to shuffle between chains Qile: the NA handling here *might* be related to the unnamed combineTCR bugs from the new rcpp con.df construction
#' @keywords internal
#' @author Ye-Lin Son Nick Borcherding
.off.the.chain <- function(dat, chain, cloneCall, check = TRUE) {
  chain1 <- toupper(chain) #to just make it easier
  if (chain1 %in% c("TRA", "TRG", "IGH")) {
    x <- 1
  } else if (chain1 %in% c("TRB", "TRD", "IGL")) {
    x <- 2
  } else {
    warning("It looks like ", chain, " does not match the available options for `chain = `")
  }
  
  if(check){
  #Adding a chain check to prevent issues with TRA + TRD/IGH data
    chain.check<- substr(str_split(dat[,"CTgene"], "_", simplify = TRUE)[,x], 1,3)
    chain.check[chain.check == "NA"] <- NA
    chain.check[chain.check == "NA;NA"] <- NA
    chain.check[chain.check == "Non"] <- NA
    any.alt.chains <- which(!is.na(chain.check) & chain.check != chain)
    if(length(any.alt.chains) > 0) {
      dat <- dat[-any.alt.chains,]
    }
  }
  
  dat[,cloneCall] <- str_split(dat[,cloneCall], "_", simplify = TRUE)[,x]
  dat[,cloneCall][dat[,cloneCall] == "NA"] <- NA
  dat[,cloneCall][dat[,cloneCall] == "NA;NA"] <- NA
 
  return(dat)
}

.padded_strings <- function(strings, max_length) {
  
  x <- lapply(strings, function(str) {
    str_len <- nchar(str)
    str <- strsplit(str, split = "")[[1]]
    if (str_len < max_length) {
      c(str, rep(NA, max_length - str_len))
    } else {
      str
    }
  })
}


#' @importFrom stringr str_split
.aa.counter <- function(input.data, cloneCall, aa.length) {
    lapply(input.data, function(x) {
        strings <- x[,cloneCall]
        strings <- do.call(c,str_split(strings, ";"))
        strings <- strings[strings != "NA"]
        strings <- na.omit(strings)
        strings <- strings[nchar(strings) < aa.length]
        strings <- .padded_strings(strings, aa.length)
        strings <- do.call(rbind, strings)
        aa.output <- apply(strings, 2, function(z) {
          summary <- as.data.frame(table(z, useNA = "always"))
        })
        res <- suppressWarnings(Reduce(function(...) merge(..., all = TRUE, by="z"), aa.output))
        colnames(res) <- c("AA", paste0("pos.", seq_len(aa.length)))
        res[seq_len(20),][is.na(res[seq_len(20),])] <- 0
        res
    }) -> res.list
  return(res.list)
}

#Pulling a color palette for visualizations
#' @importFrom grDevices hcl.colors
#' @keywords internal
.colorizer <- function(palette = "inferno", 
                        n= NULL) {
  colors <- hcl.colors(n=n, palette = palette, fixup = TRUE)
  return(colors)
}



#Remove list elements that contain all NA values
#' @keywords internal
.checkBlanks <- function(df, cloneCall) {
  count <- NULL
    for (i in seq_along(df)) {
        if (length(df[[i]][,cloneCall]) == length(which(is.na(df[[i]][,cloneCall]))) | 
            length(which(!is.na(df[[i]][,cloneCall]))) == 0 | 
            nrow(df[[i]]) == 0) {
          count <- c(i, count)
        } else {
            next()
        }
    }
  if (!is.null(count)) {
    df <- df[-count]
  }
    return(df)
}

#reshuffling df
#' @keywords internal
.groupList <- function(df, group.by) {
    df <- bind_rows(df)
    df <- split(df, df[,group.by])
    return(df)
}

# Ensure df is in list format
#' @keywords internal
.checkList <- function(df) {
  df <- tryCatch(
      {
          if (!inherits(df, "list")) { 
              df <- list(df)
          }
          df
      },
      error = function(e) {
          stop(
              "Please ensure that the input consists of at least one dataframe"
          )
      }
    )
    df
}

#' @keywords internal
.checkContigs <- function(df) {
    df <- lapply(seq_len(length(df)), function(x) {
        df[[x]] <- if(!is.data.frame(df[[x]])) as.data.frame(df[[x]]) else df[[x]]
        df[[x]][df[[x]] == ""] <- NA
        df[[x]]
    })
    df
}

#' @importFrom dplyr bind_rows
#' @keywords internal
.bound.input.return <- function(df) {
  if (is_seurat_or_se_object(df)) {
    return(.grabMeta(df))
  } 
  bind_rows(df, .id = "element.names")
}

#' @keywords internal
.list.input.return <- function(df, split.by) {
    if (is_seurat_or_se_object(df)) {
        if(is.null(split.by)){
            split.by <- "ident"
        }
        df <- .expression2List(df, split.by)
    } 
    df
}

#Get UMAP or other coordinates
#' @importFrom SingleCellExperiment reducedDim
#' @keywords internal
.get.coord <- function(sc, reduction) { 
  if (is.null(reduction)) {
    reduction <- "pca"
  }
  if (is_seurat_object(sc)) {
    coord <- sc@reductions[[reduction]]@cell.embeddings
  } else if (is_se_object(sc)) {
    coord <- reducedDim(sc, reduction)
  }
  return(coord)
}

#This is to check the single-cell expression object
#' @keywords internal
.checkSingleObject <- function(sc) {
    if (!is_seurat_or_se_object(sc)){
        stop("Object indicated is not of class 'Seurat' or 
            'SummarizedExperiment', make sure you are using
            the correct data.") 
    }
}

#This is to grab the meta data from a seurat or SCE object
#' @importFrom SingleCellExperiment colData 
#' @importFrom SeuratObject Idents
#' @importFrom methods slot
#' @keywords internal
.grabMeta <- function(sc) {
    if (is_seurat_object(sc)) {
        meta <- data.frame(sc[[]], slot(sc, "active.ident"))
        colnames(meta)[length(meta)] <- "ident"
        
    } else if (is_se_object(sc)){
        meta <- data.frame(colData(sc))
        rownames(meta) <- sc@colData@rownames
        clu <- which(colnames(meta) == "label") # as set by colLabels()
        colnames(meta)[clu] <- "ident"
    } else {
        stop("Object indicated is not of class 'Seurat' or 
            'SummarizedExperiment', make sure you are using
            the correct data.")
    }
    return(meta)
}

#This is to add the sample and ID prefixes for combineBCR()/TCR()
#' @keywords internal
.modifyBarcodes <- function(df, samples, ID) {
    out <- NULL
    for (x in seq_along(df)) {
        data <- df[[x]] 
        if (!is.null(ID)){
          data$barcode <- paste0(samples[x], "_", ID[x], "_", data$barcode)
        } else {
          data$barcode <- paste0(samples[x], "_", data$barcode)
        }
        out[[x]] <- data }
    return(out)
}

#Removing barcodes with NA recovered
#' @keywords internal
.removingNA <- function(final) {
    for(i in seq_along(final)) {
        final[[i]] <- na.omit(final[[i]])}
    return(final)
}

#Removing barcodes with > 2 clones recovered
#' @keywords internal
.removingMulti <- function(final){
    for(i in seq_along(final)) {
        final[[i]] <- filter(final[[i]], !grepl(";",CTnt))}
    return(final)
}

#Removing extra clones in barcodes with > 2 productive contigs
#' @import dplyr
#' @keywords internal
.filteringMulti <- function(x) {
    x <- x %>%
      group_by(barcode, chain) %>% 
      slice_max(n = 1, order_by = reads)
    #checking the number of productive contigs
    table <- subset(as.data.frame(table(x$barcode)), Freq > 2)
    if (nrow(table) > 0) {
      barcodes <- as.character(unique(table$Var1))
      multichain <- NULL
      for (j in seq_along(barcodes)) {
        chain <- x[x$barcode == barcodes[j],] %>% 
            group_by(barcode) %>% 
            slice_max(n = 2, order_by = reads, with_ties = FALSE)
        multichain <- rbind(multichain, chain) 
        }
    x <- subset(x, barcode %!in% barcodes)
    x <- rbind(x, multichain) 
    } 
    return(x)
}



#Filtering NA contigs out of single-cell expression object
#' @import dplyr
#' @importFrom SingleCellExperiment colData
#' @keywords internal
.filteringNA <- function(sc) {
    meta <- .grabMeta(sc)
    evalNA <- data.frame(meta[,"cloneSize"])
    colnames(evalNA) <- "indicator"
    evalNA <- evalNA %>%
        transmute(indicator = ifelse(is.na(indicator), 0, 1))
    rownames(evalNA) <- rownames(meta)
    if (inherits(x=sc, what ="cell_data_set")){
      colData(sc)[["evalNA"]]<-evalNA
      return(sc[, !is.na(sc$cloneSize)])
    }else{
      pos <- which(evalNA[,"indicator"] != 0)
      sc <- subset(sc, cells = rownames(evalNA)[pos])
      return(sc)
    }
}

#Organizing list of contigs for visualization
#' @keywords internal
.parseContigs <- function(df, i, names, cloneCall) {
    data <- df[[i]]
    data1 <- data %>% 
              group_by(data[,cloneCall]) %>%
              summarise(Abundance=n())
    colnames(data1)[1] <- cloneCall
    data1$values <- names[i]
    return(data1)
}

# This is to help sort the type of clone data to use
#' @keywords internal
.theCall <- function(df, x, check.df = TRUE) {
    x <- .convertClonecall(x)
    if(check.df) {
      if(inherits(df, "list") & !any(colnames(df[[1]]) %in% x)) {
        stop("Check the clonal variable (cloneCall) being used in the function, it does not appear in the data provided.")
      } else if (inherits(df, "data.frame") & !any(colnames(df) %in% x)) {
        stop("Check the clonal variable (cloneCall) being used in the function, it does not appear in the data provided.")
      }
    }
    return(x)
}

# helper for .theCall
.convertClonecall <- function(x) {

  clonecall_dictionary <- hash::hash(
    "gene" = "CTgene",
		"genes" = "CTgene",
		"ctgene" = "CTgene",
		"ctstrict" = "CTstrict",
		"nt" = "CTnt",
		"nucleotide" = "CTnt",
		"nucleotides" = "CTnt",
		"ctnt" = "CTnt",
		"aa" = "CTaa",
		"amino" = "CTaa",
		"ctaa" = "CTaa",
		"gene+nt" = "CTstrict",
		"strict" = "CTstrict",
		"ctstrict" = "CTstrict"
	)

	x <- tolower(x)

	if (!is.null(clonecall_dictionary[[x]])) {
		return(clonecall_dictionary[[x]])
	}
	else {
		warning("A custom variable ", x, " will be used to call clones")
		return(x)
	}
}


# Assigning positions for TCR contig data
# Used to be .parseTCR(Con.df, unique_df, data2) in v1
# but now also constructs Con.df and runs the parseTCR algorithm on it, all in Rcpp
#' @author Gloria Kraus, Nick Bormann, Nicky de Vrij, Nick Borcherding, Qile Yang
#' @keywords internal
.constructConDfAndParseTCR <- function(data2) {
  rcppConstructConDfAndParseTCR(
    data2 %>% arrange(., chain, cdr3_nt),
    unique(data2[[1]]) # 1 is the index of the barcode column
  )
}

#Assigning positions for BCR contig data
#Now assumes lambda over kappa in the context of only 2 light chains
#' @author Gloria Kraus, Nick Bormann, Nick Borcherding, Qile Yang
#' @keywords internal
.parseBCR <- function (Con.df, unique_df, data2) {
  barcodeIndex <- rcppConstructBarcodeIndex(unique_df, data2$barcode)
  for (y in seq_along(unique_df)) {
    location.i <- barcodeIndex[[y]] # *may* be wrong but should be fine. Test on old version first
    
    for (z in seq_along(location.i)) {
      where.chain <- data2[location.i[z],"chain"]
      
      if (where.chain == "IGH") {
        if(is.na(Con.df[y,"IGH"])) {
          Con.df[y,heavy_lines] <- data2[location.i[z],h_lines]
        } else {
          Con.df[y,heavy_lines] <- paste(Con.df[y, heavy_lines],
                                         data2[location.i[z],h_lines],sep=";") 
        }
      } else if (where.chain == "IGK") {
        if(is.na(Con.df[y,"IGLC"])) {
          Con.df[y,light_lines] <- data2[location.i[z],k_lines]
        } else {
          Con.df[y,light_lines] <- paste(Con.df[y, light_lines],
                                         data2[location.i[z],k_lines],sep=";") 
        }
      }else if (where.chain == "IGL") {
        if(is.na(Con.df[y,"IGLC"])) {
          Con.df[y,light_lines] <- data2[location.i[z],l_lines]
        } else {
          Con.df[y,light_lines] <- paste(Con.df[y, light_lines],
                                         data2[location.i[z],l_lines],sep=";") 
        }
      }
    }
    
  }
  return(Con.df)
}



#Producing a data frame to visualize for lengthContig()
#' @keywords internal
.lengthDF <- function(df, cloneCall, chain, group, c1, c2){
    Con.df <- NULL
    names <- names(df)
    if (identical(chain, "both")) {
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
    } else {
            for (x in seq_along(df)) {
                strings <- df[[x]][,cloneCall]
                val1 <- strings
                for (i in seq_along(val1)) {
                    if (grepl(";", val1[i]) == TRUE) {
                        val1[i] <- str_split(val1[i], ";", simplify = TRUE)[1] 
                    } else { next() } }
                chain1 <- nchar(val1)
                if (!is.null(group)) {
                    cols1 <- df[[x]][,group]
                    data1 <- data.frame(chain1, val1, names[x], c1, cols1)
                    colnames(data1)<-c("length","CT","values","chain",group)
                }else if (is.null(group)){
                    data1 <- data.frame(chain1, val1, names[x], c1)
                    colnames(data1) <- c("length", "CT", "values", "chain")
                data <- na.omit(data1)
                data <- subset(data, CT != "NA" & CT != "")
                Con.df<- rbind.data.frame(Con.df, data) }}
    }
    return(Con.df)
}

#General combination of nucleotide, aa, and gene sequences for T/B cells
#' @keywords internal
.assignCT <- function(cellType, Con.df) {
    if (cellType == "T") {
        Con.df$CTgene <- paste(Con.df$TCR1, Con.df$TCR2, sep="_")
        Con.df$CTnt <- paste(Con.df$cdr3_nt1, Con.df$cdr3_nt2, sep="_")
        Con.df$CTaa <- paste(Con.df$cdr3_aa1, Con.df$cdr3_aa2, sep="_")
        Con.df$CTstrict <- paste0(Con.df$TCR1, ";", Con.df$cdr3_nt1, "_",
            Con.df$TCR2, ";", Con.df$cdr3_nt2)
    } else { # assume cellType = B
        Con.df$CTgene <- paste(Con.df$IGH, Con.df$IGLC, sep="_")
        Con.df$CTnt <- paste(Con.df$cdr3_nt1, Con.df$cdr3_nt2, sep="_")
        Con.df$CTaa <- paste(Con.df$cdr3_aa1, Con.df$cdr3_aa2, sep="_")
    }
    return(Con.df)
}


#Sorting the V/D/J/C gene sequences for T and B cells
#' @importFrom stringr str_c str_replace_na
#' @importFrom dplyr bind_rows
#' @keywords internal
.makeGenes <- function(cellType, data2, chain1, chain2) {
    if(cellType %in% c("T")) {
        data2 <- data2 %>% 
            mutate(TCR1 = ifelse(chain %in% c("TRA", "TRG"), 
                  str_c(str_replace_na(v_gene),  str_replace_na(j_gene), str_replace_na(c_gene), sep = "."), NA)) %>%
            mutate(TCR2 = ifelse(chain %in% c("TRB", "TRD"), 
                  str_c(str_replace_na(v_gene), str_replace_na(d_gene),  str_replace_na(j_gene),  str_replace_na(c_gene), sep = "."), NA))
    } else if (cellType %in% c("B")) {
        heavy <- data2[data2$chain == "IGH",] %>% 
          mutate(IGHct = str_c(str_replace_na(v_gene), str_replace_na(d_gene),  str_replace_na(j_gene),  str_replace_na(c_gene), sep = "."))
        kappa <- data2[data2$chain == "IGK",] %>% 
          mutate(IGKct = str_c(str_replace_na(v_gene),  str_replace_na(j_gene), str_replace_na(c_gene), sep = ".")) 
        lambda <- data2[data2$chain == "IGL",] %>%
          mutate(IGLct = str_c(str_replace_na(v_gene),  str_replace_na(j_gene), str_replace_na(c_gene), sep = "."))
        data2 <- bind_rows(heavy, kappa, lambda)
    }
    data2
}

#Getting the minimum repertoire size for diversity boots
#' @keywords internal
.short.check <- function(df, cloneCall) {
  min <- c()
  for (x in seq_along(df)) {
    min.tmp <- length(which(!is.na(unique(df[[x]][,cloneCall]))))
    min <- c(min.tmp, min)
  }
  min <- min(min)
  return(min)
}

# check if object is a dataframe or list of dataframes
#' @keywords internal
is_df_or_list_of_df <- function(x) {
    if (is.data.frame(x)) {
        return(TRUE)
    } else if (is.list(x)) {
        if (length(x) == 0) {
            return(FALSE)
        }
        return(all(sapply(x, is.data.frame)))
    } else {
        return(FALSE)
    }
}

#Pulling meta data 
#' @keywords internal
.expression2List <- function(sc, split.by) {
  if (!inherits(x=sc, what ="Seurat") & 
      !inherits(x=sc, what ="SummarizedExperiment")) {
    stop("Use a Seurat or SCE object to convert into a list")
  }
  meta <- .grabMeta(sc)
  if(is.null(split.by)){
    split.by <- "ident"
  }
  unique <- str_sort(as.character(unique(meta[,split.by])), numeric = TRUE)
  df <- NULL
  for (i in seq_along(unique)) {
    subset <- subset(meta, meta[,split.by] == unique[i])
    subset <- subset(subset, !is.na(cloneSize))
    df[[i]] <- subset
  }
  names(df) <- unique
  return(df)
}

#Making lists for single-cell object, check blanks and apply chain filter
#' @keywords internal
.data.wrangle <- function(df, split.by, cloneCall, chain) {
  df <- .list.input.return(df, split.by)
  df <- .checkBlanks(df, cloneCall)
  for (i in seq_along(df)) {
    if (chain != "both") {
      df[[i]] <- .off.the.chain(df[[i]], chain, cloneCall)
    }
  }
  return(df)
}

# Calculates the normalized Levenshtein Distance between the contig 
# nucleotide sequence.
#' @importFrom stringdist stringdist
#' @importFrom igraph graph_from_data_frame components graph_from_edgelist
#' @importFrom dplyr bind_rows
#' @keywords internal
.lvCompare <- function(dictionary, gene, chain, threshold, exportGraph = FALSE) {
  overlap <- NULL
  out <- NULL
  dictionary[dictionary == "None"] <- NA
  dictionary$v.gene <- stringr::str_split(dictionary[,gene], "[.]", simplify = TRUE)[,1]
  tmp <- na.omit(unique(dictionary[,c(chain, "v.gene")]))
  #chunking the sequences for distance by v.gene
  unique.v <- na.omit(unique(tmp$v.gene))
  
  edge.list <- lapply(unique.v, function(y) {
  #for(y in unique.v) {
    secondary.list <- list()
    filtered_df <- tmp %>% filter(v.gene == y)
    filtered_df <- filtered_df[!is.na(filtered_df[,chain]),]
    nucleotides <- unique(filtered_df[,chain])
    if (length(nucleotides) > 1) {
      chain_col_number <- 1
      nucleotide_lengths <- nchar(nucleotides)
      # Pre-allocate list
      list <- vector("list", length = length(nucleotides))
      
      for (i in seq_len((length(nucleotides) - 1))) {
        temp_list <- vector("list", length = length(nucleotides) - i)
        
        idx_i <- tmp[, chain_col_number] == nucleotides[i] & tmp[,2] == y
        len_i <- nucleotide_lengths[i]
        
        for (j in (i + 1):length(nucleotides)) {
          distance <- stringdist::stringdist(nucleotides[i], nucleotides[j], method = "lv")
          distance_norm <- 1 - distance / ((len_i + nucleotide_lengths[j]) / 2)
          
          if (!is.na(distance_norm) & distance_norm >= threshold) {
            idx_j <- tmp[, chain_col_number] == nucleotides[j]
            stored_positions <- idx_j & !idx_i
            
            if(any(stored_positions)) {
              # Store this pair in the edge list that is not the same chain
              temp_list[[j - i]] <- data.frame(from = which(idx_i), 
                                               to = which(stored_positions))
            }
          }
        }
        list[[i]] <- do.call(rbind, temp_list)  # Collapsing all data.frames in temp_list
      }
      
      #Remove any NULL or 0 list elements
      if(length(list) > 0) {
        list <-  list[-which(unlist(lapply(list, is.null)))]
        list <-  list[lapply(list,length)>0]
        list <- bind_rows(list) %>% as.data.frame()
        secondary.list[[i]] <- list
      }
      #Remove any NULL or 0 list elements
      if(length(secondary.list) > 0) {
        secondary.list <-  secondary.list[-which(unlist(lapply(secondary.list, is.null)))]
        secondary.list <-  secondary.list[lapply(secondary.list,length)>0]
      }
    }
    secondary.list
  })
  edge.list <- bind_rows(edge.list)
  if(exportGraph) {
    edge.list[,1] <- tmp[,1][edge.list[,1]]
    edge.list[,2] <- tmp[,1][edge.list[,2]]
    graph <- graph_from_edgelist(as.matrix(edge.list))
    return(graph)
  }
  
  if (nrow(edge.list) > 0) { 
    edge.list <- unique(edge.list)
    g <- graph_from_data_frame(edge.list)
    components <- igraph::components(g, mode = c("weak"))
    out <- data.frame("cluster" = components$membership, 
                      "filtered" = names(components$membership))
    filter <- which(table(out$cluster) > 1)
    out <- subset(out, cluster %in% filter)
    if(nrow(out) > 1) {
      out$cluster <- paste0(gene, ":Cluster", ".", out$cluster)
      out$filtered <- tmp[,1][as.numeric(out$filtered)]
      uni_IG <- as.data.frame(unique(tmp[,1][tmp[,1] %!in% out$filtered]))
    }
  } else {
    uni_IG <- as.data.frame(unique(tmp[,1]))
  }
  colnames(uni_IG)[1] <- "filtered"
  if (nrow(uni_IG) > 0) {
    uni_IG$cluster <- paste0(gene, ".", seq_len(nrow(uni_IG)))
  }
  output <- rbind.data.frame(out, uni_IG)
  colnames(output) <- c("Cluster", "clone")
  return(output)
}

