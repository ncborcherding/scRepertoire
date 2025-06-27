#-------------------------------------------
#-------------Background Functions-----------
#-------------------------------------------
# utility functions use camelCase 

.toCapitilize <- function(name) {
  gsub("([\\w])([\\w]+)", "\\U\\1\\L\\2", name, perl = TRUE)
}

.orderingFunction <- function(vector,
                               group.by,
                               data.frame) {
  if (length(vector) == 1 && vector == "alphanumeric") {
    unique_elements <- as.character(unique(data.frame[, group.by]))
    sorted_levels <- unique_elements[
      order(gsub("[0-9]", "", unique_elements),
            as.numeric(gsub("[^0-9]", "", unique_elements)))
    ]
    # Apply the sorted levels to the factor
    data.frame[, group.by] <- factor(data.frame[, group.by], levels = sorted_levels)
  } else {
    data.frame[, group.by] <- factor(data.frame[, group.by], levels = vector)
  }
  
  return(data.frame)
}

# Get position of chain
.chainPositionParser <- function(chain) {
  chain1 <- toupper(chain) #to just make it easier
  if (chain1 %in% c("TRA", "TRG", "IGH")) {
    x <- 1
  } else if (chain1 %in% c("TRB", "TRD", "IGL", "IGK")) {
    x <- 2
  } else {
    # It's good practice to use stop() for fatal errors
    stop("'", chain, "' is not a valid entry for the `chain` argument.
         Please use 'TRA', 'TRG', 'IGH', 'TRB', 'TRD', 'IGL', 'IGK'.")
  }
  return(x)
}

#Use to shuffle between chains
#' @keywords internal
#' @author Ye-Lin Son, Nick Borcherding
.offTheChain <- function(dat, chain, cloneCall, check = TRUE) {
  # Get position of chain
  x <- .chainPositionParser(chain)
  
  if (check) {
    split_genes <- strsplit(dat[, "CTgene"], "_")
    gene_elements <- sapply(split_genes, `[`, x)
    
    chain.check <- substr(gene_elements, 1, 3)
    
    # The following lines for handling NAs remain the same
    chain.check[chain.check == "NA"] <- NA
    chain.check[is.na(chain.check)] <- NA #
    chain.check[chain.check == "Non"] <- NA
    
    any.alt.chains <- which(!is.na(chain.check) & chain.check != chain)
    
    if (length(any.alt.chains) > 0) {
      dat <- dat[-any.alt.chains, ]
      split_clones <- strsplit(dat[, cloneCall], "_")
    } else {
      split_clones <- strsplit(dat[, cloneCall], "_")
    }
  } else {
    split_clones <- strsplit(dat[, cloneCall], "_")
  }
  
  # Extract the relevant part of the cloneCall
  dat[, cloneCall] <- sapply(split_clones, `[`, x)
  
  # Handle NA-like values.
  dat[, cloneCall][dat[, cloneCall] == "NA"] <- NA
  dat[, cloneCall][is.na(dat[, cloneCall])] <- NA
  
  return(dat)
}

.cloneCounter <- function(meta,
                           group.by,
                           cloneCall) {
  
  meta_filtered <- na.omit(meta[, c(group.by, cloneCall)])
  clone.table <- as.data.frame(table(meta_filtered), responseName = "n")
  clone.table <- clone.table[order(clone.table$n, decreasing = TRUE), ]
  clone.table$group.sum <- ave(clone.table$n, clone.table[, group.by], FUN = sum)
  clone.table$clone.sum <- ave(clone.table$n, clone.table[, cloneCall], FUN = sum)
  rownames(clone.table) <- NULL
  return(clone.table)
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
    # First, check if there are no rows
    if (nrow(df[[i]]) == 0) {
      count <- c(i, count)
    } 
    # If there are rows, then proceed with blank checks
    else if (length(df[[i]][,cloneCall]) == length(which(is.na(df[[i]][,cloneCall]))) |
             length(which(!is.na(df[[i]][,cloneCall]))) == 0) {
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
    df <- do.call(rbind, df)
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
              "Please ensure that the input consists of at least one data frame"
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

#' @keywords internal
.bound.input.return <- function(df) {
  if (.is.seurat.or.se.object(df)) {
    return(.grabMeta(df))
  } 
  bind_rows(df, .id = "element.names")
}

#' @keywords internal
.list.input.return <- function(df, split.by) {
    if (.is.seurat.or.se.object(df)) {
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
.getCoord <- function(sc, reduction) { 
  if (is.null(reduction)) {
    reduction <- "pca"
  }
  if (.is.seurat.object(sc)) {
    if (!reduction %in% names(sc@reductions)) {
      stop(paste("Reduction", reduction, "not found in Seurat object."))
    }
    coord <- sc@reductions[[reduction]]@cell.embeddings
  } else if (.is.se.object(sc)) {
    if (!reduction %in% reducedDimNames(sc)) {
      stop(paste("Reduction", reduction, "not found in SingleCellExperiment object."))
    }
    coord <- reducedDim(sc, reduction)
  }
  return(coord)
}

#This is to check the single-cell expression object
#' @keywords internal
.checkSingleObject <- function(sc) {
    if (!.is.seurat.or.se.object(sc)){
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
  messString <- c("Meta data contains an 'ident' column and will likely result
                  in errors downstream.")
    if (.is.seurat.object(sc)) {
        meta <- data.frame(sc[[]], slot(sc, "active.ident"))
        if("ident" %in% colnames(meta)) {
          message(messString)
        }
        colnames(meta)[length(meta)] <- "ident"
    } else if (.is.se.object(sc)){
        meta <- data.frame(colData(sc))
        if("ident" %in% colnames(meta)) {
          message(messString)
        }
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
#' @keywords internal
.filteringMulti <- function(x) {
  x <- as.data.frame(x)
  x <- x[order(x$reads, decreasing = TRUE), ]
  x <- x[!duplicated(x[, c("barcode", "chain")]), ]
  contig_counts <- as.data.frame(table(x$barcode))
  table <- subset(contig_counts, Freq > 2)
  if (nrow(table) > 0) {
    barcodes <- as.character(unique(table$Var1))
    multichain <- NULL
    multi_subset <- subset(x, barcode %in% barcodes)
    multi_subset <- multi_subset[order(multi_subset$barcode, multi_subset$reads, decreasing = TRUE), ]
    multi_subset$rank <- ave(multi_subset$reads, multi_subset$barcode, FUN = seq_along)
    multichain <- subset(multi_subset, rank <= 2)
    multichain$rank <- NULL
    x <- subset(x, !barcode %in% barcodes)
    x <- rbind(x, multichain)
  }
  x <- x[order(x$barcode, x$chain), ]
  rownames(x) <- NULL 
  return(x)
}


#Filtering NA contigs out of single-cell expression object
#' @importFrom SingleCellExperiment colData
#' @keywords internal
.filteringNA <- function(sc) {
  meta <- .grabMeta(sc)
  evalNA <- data.frame(indicator = meta[, "cloneSize"], 
                       row.names = rownames(meta))
  evalNA$indicator <- ifelse(is.na(evalNA$indicator), 0, 1)
  if (.is.se.object(sc)) {
    colData(sc)[["evalNA"]] <- evalNA
    return(sc[, !is.na(sc$cloneSize)])
  } else {
    pos <- which(evalNA[, "indicator"] != 0)
    sc <- subset(sc, cells = rownames(evalNA)[pos])
    return(sc)
  }
}

#Organizing list of contigs for visualization
#' @keywords internal
#' @author Justin Reimertz
.parseContigs <- function(df, i, names, cloneCall) {
	data <- df[[i]] %>% 
		# Count the abundance of each clone per specified cloneCall category and
		# cell type
		group_by("{cloneCall}" := df[[i]][[cloneCall]]) %>%
		# Add sample ids to df
		mutate(values = names[i], Abundance = n()) %>%
		select(all_of(cloneCall), values, Abundance) %>%
		ungroup() 
	return(data)
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

# helper for .theCall # Qile: on second thought - converting to x to lowercase may be a bad idea...
.convertClonecall <- function(x) {

  clonecall_dictionary <- list(
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

	if (!is.null(clonecall_dictionary[[tolower(x)]])) {
		return(clonecall_dictionary[[tolower(x)]])
	} else {
		message("A custom variable ", x, " will be used to call clones")
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
    dplyr::arrange(data2, chain, cdr3_nt),
    uniqueData2Barcodes = unique(data2$barcode)
  )
}

# Assigning positions for BCR contig data
# Now assumes lambda over kappa in the context of only 2 light chains
#' @author Gloria Kraus, Nick Bormann, Nick Borcherding, Qile Yang
#' @keywords internal
.parseBCR <- function (Con.df, unique_df, data2) {
  barcodeIndex <- rcppConstructBarcodeIndex(unique_df, data2$barcode)
  for (y in seq_along(unique_df)) {
    location.i <- barcodeIndex[[y]]
    
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



# Producing a data frame to visualize for lengthContig()
#' @keywords internal
.lengthDF <- function(df, cloneCall, chain, group) {
  Con.df <- NULL
  names <- names(df)
  
  if (identical(chain, "both")) {
    for (i in seq_along(df)) {
      clones <- gsub("_NA", "", df[[i]][, cloneCall])
      clones <- gsub("NA_", "", clones)
      length <- nchar(gsub("_", "", clones))
      val <- df[[i]][, cloneCall]
      
      if (!is.null(group)) {
        cols <- df[[i]][, group]
        data <- na.omit(data.frame(length, val, cols, names[i]))
        colnames(data) <- c("length", "CT", group, "values")
        Con.df <- rbind(Con.df, data) 
      } else {
        data <- na.omit(data.frame(length, val, names[i]))
        colnames(data) <- c("length", "CT", "values")
        Con.df <- rbind(Con.df, data)
      }
    }
  } else {
    for (x in seq_along(df)) {
      strings <- df[[x]][, cloneCall]
      val1 <- strings
      split_strings <- strsplit(val1, ";")
      val1 <- sapply(split_strings, `[`, 1)
      chain1 <- nchar(val1)
      
      if (!is.null(group)) {
        cols1 <- df[[x]][, group]
        data1 <- data.frame(chain1, val1, names[x], chain, cols1)
        colnames(data1) <- c("length", "CT", "values", "chain", group)
      } else { 
        data1 <- data.frame(chain1, val1, names[x], chain)
        colnames(data1) <- c("length", "CT", "values", "chain")
      }
      
      data <- na.omit(data1)
      data <- data[data$CT != "NA" & data$CT != "" & !is.na(data$CT), ]
      Con.df <- rbind(Con.df, data)
    }
  }
  return(Con.df)
}
# General combination of nucleotide, aa, and gene sequences for T/B cells
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


# Sorting the V/D/J/C gene sequences for T and B cells
#' @keywords internal
.makeGenes <- function(cellType, data2, chain1, chain2) {

  replace_na_with_empty <- function(x) {
    x[is.na(x)] <- ""
    return(x)
  }
  
  if (cellType == "T") {
    # Filter for T-cell chains
    data2 <- data2[data2$chain %in% c("TRA", "TRB", "TRG", "TRD"), ]
    
    # Create TCR1 column for TRA/TRG chains
    is_tcr1_chain <- data2$chain %in% c("TRA", "TRG")
    data2$TCR1 <- ifelse(
      is_tcr1_chain,
      paste(
        replace_na_with_empty(data2$v_gene),
        replace_na_with_empty(data2$j_gene),
        replace_na_with_empty(data2$c_gene),
        sep = "."
      ),
      NA
    )
    
    # Create TCR2 column for TRB/TRD chains
    is_tcr2_chain <- data2$chain %in% c("TRB", "TRD")
    data2$TCR2 <- ifelse(
      is_tcr2_chain,
      paste(
        replace_na_with_empty(data2$v_gene),
        replace_na_with_empty(data2$d_gene),
        replace_na_with_empty(data2$j_gene),
        replace_na_with_empty(data2$c_gene),
        sep = "."
      ),
      NA
    )
    
  } else if (cellType == "B") {
    # Filter for B-cell chains
    data2 <- data2[data2$chain %in% c("IGH", "IGK", "IGL"), ]
    
    data2$IGHct <- NA
    data2$IGKct <- NA
    data2$IGLct <- NA
    
    # Populate IGHct for IGH chains
    is_heavy <- data2$chain == "IGH" & !is.na(data2$chain)
    data2$IGHct[is_heavy] <- paste(
      replace_na_with_empty(data2$v_gene[is_heavy]),
      replace_na_with_empty(data2$d_gene[is_heavy]),
      replace_na_with_empty(data2$j_gene[is_heavy]),
      replace_na_with_empty(data2$c_gene[is_heavy]),
      sep = "."
    )
    
    # Populate IGKct for IGK chains
    is_kappa <- data2$chain == "IGK" & !is.na(data2$chain)
    data2$IGKct[is_kappa] <- paste(
      replace_na_with_empty(data2$v_gene[is_kappa]),
      replace_na_with_empty(data2$j_gene[is_kappa]),
      replace_na_with_empty(data2$c_gene[is_kappa]),
      sep = "."
    )
    
    # Populate IGLct for IGL chains
    is_lambda <- data2$chain == "IGL" & !is.na(data2$chain)
    data2$IGLct[is_lambda] <- paste(
      replace_na_with_empty(data2$v_gene[is_lambda]),
      replace_na_with_empty(data2$j_gene[is_lambda]),
      replace_na_with_empty(data2$c_gene[is_lambda]),
      sep = "."
    )
  }
  return(data2)
}




# Pulling meta data 
#' @keywords internal
.expression2List <- function(sc, split.by) {
  .checkSingleObject(sc)
  meta <- .grabMeta(sc)
  if(is.null(split.by)){
    split.by <- "ident"
  }
  unique_elements <- as.character(unique(meta[, split.by]))
  sorted_unique <- unique_elements[
    order(gsub("[0-9]", "", unique_elements),
          as.numeric(gsub("[^0-9]", "", unique_elements)))
  ]
  df <- list()
  
  for (i in seq_along(sorted_unique)) {
    current_level <- sorted_unique[i]
    
    # Filter metadata for the current level
    temp_df <- meta[meta[, split.by] == current_level & !is.na(meta[, split.by]), ]
    
    # Filter out rows where cloneSize is NA, checking if the column exists first
    if ("cloneSize" %in% colnames(temp_df)) {
      temp_df <- temp_df[!is.na(temp_df$cloneSize), ]
    }
    
    df[[i]] <- temp_df
  }
  
  names(df) <- sorted_unique
  return(df)
}

# Making lists for single-cell object, check blanks and apply chain filter
#' @keywords internal
.dataWrangle <- function(df, split.by, cloneCall, chain) {
  df <- .list.input.return(df, split.by)
  df <- .checkBlanks(df, cloneCall)
  for (i in seq_along(df)) {
    if (chain != "both") {
      df[[i]] <- .offTheChain(df[[i]], chain, cloneCall)
    }
  }
  return(df)
}


#-------------------------------------------
#-----------------Type Checks---------------
#-------------------------------------------
# Separate from utility functions with dot designations
.is.seurat.object <- function(obj) inherits(obj, "Seurat")

.is.se.object <- function(obj) inherits(obj, "SummarizedExperiment")

.is.seurat.or.se.object <- function(obj) {
  .is.seurat.object(obj) || .is.se.object(obj)
}

.is.named.numeric <- function(obj) {
  is.numeric(obj) && !is.null(names(obj))
}

.is.df.or.list.of.df <- function(x) {
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



