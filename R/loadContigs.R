#' Loading the contigs derived from single-cell sequencing
#'
#' This function generates a contig list and formats the data to allow for 
#' function with  \code{\link{combineTCR}} or \code{\link{combineBCR}}. If 
#' using data derived from filtered outputs of 10X Genomics, there is no 
#' need to use this function as the data is already compatible. The function
#' assumes if listing multiple directories, there are distinct outputs with
#' unmodified file names in them. The files that this function parses includes: 
#' WAT3R = "barcode_results.csv", 10X =  "filtered_contig_annotation.csv", 
#' AIRR = "airr_rearrangement.tsv", TRUST4 = "barcode_report.tsv", and 
#' BD = "Contigs_AIRR.tsv".
#' 
#' @examples
#' \dontrun{
#' dir <- c("./Sample1/outs/")
#' contig.list <- loadContigs(dir, format = "10X")
#' }
#' 
#' TRUST4 <- read.csv("https://www.borch.dev/uploads/contigs/TRUST4_contigs.csv")
#' contig.list <- loadContigs(TRUST4, format = "TRUST4")
#' 
#' BD <- read.csv("https://www.borch.dev/uploads/contigs/BD_contigs.csv")
#' contig.list <- loadContigs(BD, format = "BD")
#' 
#' WAT3R <- read.csv("https://www.borch.dev/uploads/contigs/WAT3R_contigs.csv")
#' contig.list <- loadContigs(WAT3R, format = "WAT3R")
#' 
#' @param dir The directory in which contigs are located or a list with contig elements
#' @param format The format of the single-cell contig, currently supporting: 
#' "10X", "AIRR", "TRUST4", "BD", and "WAT3R"
#' @export
#' @return List of contigs for further processing in scRepertoire
loadContigs <- function(dir, 
                        format = "10X") {
  #Loading from directory, recursively
  if (inherits(x=dir, what ="character")) {
    format.list <- list("WAT3R" = "barcode_results.csv", 
                        "10X" =  "filtered_contig_annotation.csv", 
                        "AIRR" = "airr_rearrangement.tsv", 
                        "TRUST4" = "barcode_report.tsv", 
                        "BD" = "Contigs_AIRR.tsv")
        file.pattern <- format.list[[format]]
        contig.files <- list.files(dir, file.pattern, recursive = TRUE, full.names = TRUE)
        if (format %in% c("10X", "WAT3R")) {
          df <- lapply(contig.files, read.csv) 
        } else if (format %!in% c("10X", "WAT3R")) {
          df <- lapply(contig.files, read.delim)
        }
  #Already loaded list or data frame
  } else if (inherits(x=dir, what ="list") | inherits(x=dir, what ="data.frame")) {
    df <- checkList(dir)
  }
        
        if (format %in% c("10X", "WAT3R")) {
            if (format =="WAT3R") {
                df <- parseWAT3R(df)
            } else {
                df <- parse10x(df)
            }
        } else {
            if (format =="TRUST4") {
                df <- parseTRUST4(df) 
            } else if (format == "AIRR") {
                df <- parseAIRR(df)
            } else {
              df <- parseBD(df)
            }
        }
    return(df)
}

#Formats TRUST4 data
#' @importFrom stringr str_split
parseTRUST4 <- function(df) {
    for (i in seq_along(df)) {
        colnames(df[[i]])[1] <- "barcode"
        df[[i]][df[[i]] == "*"] <- NA
        
        if(length(which(is.na(df[[i]]$chain1))) == length(df[[i]]$chain1)) {
          chain2 <- matrix(ncol = 7, nrow = length(df[[i]]$chain1))
        } else {
          chain2 <- str_split(df[[i]]$chain1, ",", simplify = TRUE)[,seq_len(7)]
          chain2[chain2 == "*"] <- "None"
        }
        colnames(chain2) <- c("v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "reads")
        chain2 <- data.frame(barcode = df[[i]][,1], chain2)
        
        if(length(which(is.na(df[[i]]$chain2))) == length(df[[i]]$chain2)) {
          chain1 <- matrix(ncol = 7, nrow = length(df[[i]]$chain2))
        } else {
          chain1 <- str_split(df[[i]]$chain2, ",", simplify = TRUE)[,seq_len(7)]
          chain1[chain1 == "*"] <- "None"
        }
        colnames(chain1) <- c("v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "reads")
        chain1 <- data.frame(barcode = df[[i]][,1], chain1)
        data2 <- rbind(chain1, chain2)
        data2[data2 == ""] <- NA
        df[[i]] <- data2
    }
    df <- chain.parser(df)
    return(df)
}

#Formats wat3r data
#' @author Kyle Romine, Nick Borcherding
parseWAT3R <- function(df) {
    for (i in seq_along(df)) {
        colnames(df[[i]])[1] <- "barcode"
        df[[i]][df[[i]] == ""] <- NA
        chain2 <- df[[i]][,c(1,8,9,10,4,3,7,5)]
        chain2 <- data.frame(chain2[,1], chain = "TRB", chain2[,2:4], c_gene = NA, chain2[,5:8])
        colnames(chain2) <- c("barcode", "chain", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "reads", "umis")
        
        #TRA Chain 1
        chain1 <-  df[[i]][,c(1,16,17,12,11,15,13)]
        chain1 <- data.frame(chain1[,1], chain = "TRA",chain1[,2], d_gene = NA, chain1[,3], c_gene = NA, chain1[,4:7])
        colnames(chain1) <- c("barcode", "chain", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "reads", "umis")
        data2 <- rbind(chain1, chain2)
        data2[data2 == ""] <- NA
        
        #TRA Chain 2
        chain3 <-  df[[i]][,c(1,23,24,19,18,22,20)]
        chain3 <- data.frame(chain3[,1], chain = "TRA",chain3[,2],  d_gene = NA, chain3[,3], c_gene = NA, chain3[,4:7])
        colnames(chain3) <- c("barcode", "chain", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "reads", "umis")
        data2 <- rbind(chain1, chain2, chain3)
        data2[data2 == ""] <- NA
        df[[i]] <- data2
        df[[i]] <- df[[i]][with(df[[i]], order(reads, chain)),]
    }
    return(df)
}

#Formats AIRR data
parseAIRR <- function(df) {
    for (i in seq_along(df)) {
        df[[i]] <- df[[i]][,c("cell_id", "locus", "consensus_count", "v_call", "d_call", "j_call", "c_call", "junction", "junction_aa")]
        colnames(df[[i]]) <- c("barcode", "chain", "reads", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3")
        df[[i]] <- df[[i]][with(df[[i]], order(reads, chain)),]
    }
    return(df)
}

#Loads 10x data
parse10x <- function(df) {
    for (i in seq_along(df)) {
        df[[i]] <- subset(df[[i]], chain != "Multi")
        df[[i]] <- subset(df[[i]], productive %in% c(TRUE, "TRUE", "True", "true"))
        if (nrow(df[[i]]) == 0) { stop(
            "There are 0 contigs after internal filtering -
            check the contig list to see if any issues exist 
            for productive chains", call. = FALSE) }
        df[[i]] <- subset(df[[i]], cdr3 != "None")
        df[[i]][df[[i]] == ""] <- NA
        df[[i]] <- df[[i]][with(df[[i]], order(reads, chain)),]
    }
    return(df)
}
#Loads BD AIRR
parseBD <- function(df) {
  for (i in seq_along(df)) {
    df[[i]] <- df[[i]][,c(1,2,20,25,30, 35, 48,49,4)]
    colnames(df[[i]]) <- c("barcode", "chain", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "reads")
    df[[i]] <- df[[i]][with(df[[i]], order(reads, chain)),]
  }
  return(df)
}
#Grabs the chain info from v_gene
chain.parser <- function(df) {
  for (i in seq_along(df)) {
    df[[i]]$chain <- substr(df[[i]][,"v_gene"],1,3)
  }
  return(df)
}
    
