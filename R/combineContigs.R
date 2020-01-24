#df indicates the list of data frames of the contigs
#samples either single or multiple samples
#ID are sub designated for labeling of samples
#cor can be used to produce a string of sample numbers to correspond to the df list, e.g. df <- list(Sample1_Peripheral, Sample2_Tumor, Sample1_Peripheral) -- cor = c(1,2,1)
#' @export
combineContigs <- function(df,
                           samples = NULL,
                           ID = NULL, cells = c("T-AB", "T-GD", "B")) {
    df <- if(class(df) != "list") list(df) else df
    require(dplyr)
    out <- NULL
    final <- NULL
    count <- length(unlist(strsplit(df[[1]]$barcode[1], "[-]")))
    count2 <- length(unlist(strsplit(df[[1]]$barcode[1], "[_]")))
    if (count > 2 | count2 > 2) {
        stop("Seems to be an error in the naming of the contigs, ensure the barcodes are labeled like, AAACGGGAGATGGCGT-1 or AAACGGGAGATGGCGT, use stripBarcode to get the basic format", call.=F)
    } else if (length(df) != length(samples) | length(df) != length(ID)) {
        stop("Make sure the sample and ID labels match the length of the list of data frames (df).", call. = F)
    } else {
        if (cells == "T-AB" | cells == "T-GD") {
            if (cells == "T-AB") {
                chain1 <- "TRA"
                chain2 <- "TRB"
            }
            else if (cells == "T-GD") {
                chain1 <- "TRG"
                chain2 <- "TRD"
            }
            for (i in seq_along(df)) {
                df[[i]] <- subset(df[[i]], chain != "Multi")
                df[[i]] <- subset(df[[i]], chain == chain1 | chain == chain2)
                df[[i]] <- subset(df[[i]], productive == T | productive == "TRUE")
                df[[i]]$sample <- samples[i]
                df[[i]]$ID <- ID[i]
                if (nrow(df[[i]]) == 0) {
                    stop("Check some hypotenuses, Captain. There are 0 contigs after filtering for celltype.", call. = F)
                }
            }
        } else if(cells == "B") {
            for (i in seq_along(df)) {
                df[[i]] <- subset(df[[i]], chain != "Multi")
                df[[i]] <- subset(df[[i]], chain == "IGH" | chain == "IGK" | chain == "IGL")
                df[[i]] <- subset(df[[i]], productive == T | productive == "TRUE")
            }
        }

        for (x in seq_along(df)) {
            data <- df[[x]]
            data$barcode <- paste(samples[x], "_", ID[x], "_", data$barcode, sep="")
            out[[x]] <- data

        }
    }

    if (cells == "T-AB" | cells ==  "T-GD") {
        for (i in seq_along(out)) {

            data2 <- out[[i]]
            data2 <- data2 %>%
                mutate(TCR1 = ifelse(chain == chain1, paste(with(data2, interaction(v_gene,  j_gene, c_gene))), NA)) %>%
                mutate(TCR2 = ifelse(chain == chain2, paste(with(data2, interaction(v_gene,  j_gene, d_gene, c_gene))), NA))
            unique_df <- unique(data2$barcode)
            Con.df <- data.frame(matrix(NA, length(unique_df), 7))

            colnames(Con.df) <- c("barcode","TCR1", "cdr3_aa1", "cdr3_nt1", "TCR2", "cdr3_aa2", "cdr3_nt2")
            Con.df$barcode <- unique_df
            y <- NULL
            for (y in 1:length(unique_df)){
                barcode.i <- Con.df$barcode[y]
                location.i <- which(barcode.i == data2$barcode)
                if (length(location.i) == 2){
                    if (is.na(data2[location.i[1],c("TCR1")])) {
                        Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")] <- data2[location.i[1],c("TCR2", "cdr3", "cdr3_nt")]
                        Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- data2[location.i[2],c("TCR1", "cdr3", "cdr3_nt")]
                    } else {
                        Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- data2[location.i[1],c("TCR1", "cdr3", "cdr3_nt")]
                        Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")] <- data2[location.i[2],c("TCR2", "cdr3", "cdr3_nt")]
                    }
                } else if (length(location.i) == 3) { # if there are more than 2 chains for each cell
          if (is.na(data2[location.i[1],c("TCR1")])) { # if the 1st location is occupied by TRB
            Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")] <- data2[location.i[1],c("TCR2", "cdr3", "cdr3_nt")] # add TRB info to the TCR2 related columns
            if (is.na(data2[location.i[2],c("TCR1")])) { # if the 2nd location is occupied by TRB 
              TRdf <- paste(Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")],data2[location.i[2],c("TCR2", "cdr3", "cdr3_nt")],sep=";") # paste TRBs from 1st and 2nd locations
              Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")] <- TRdf # add the combination of 2 TRB chains to the TCR2 slot
              Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- data2[location.i[3],c("TCR1", "cdr3", "cdr3_nt")] # add TRA from the 3rd location to the TCR1 related columns 
            } else { # if the 2nd location is occupied by TRA
              Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- data2[location.i[1],c("TCR1", "cdr3", "cdr3_nt")] # add TRA info to the TCR1 related columns 
              if (is.na(data2[location.i[3],c("TCR1")])) { # if the 3rd location is occupied by TRB
                TRdf <- paste(Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")],data2[location.i[3],c("TCR2", "cdr3", "cdr3_nt")],sep=";") # paste TRB from 1st and 3rd locations
                Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")] <- TRdf # add combined TRB relatedn info to the TCR2 related columns
              } else { # if the 3rd location is occupied by TRA
                TRdf <- paste(Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")],data2[location.i[3],c("TCR1", "cdr3", "cdr3_nt")],sep=";") # paste TRAs from 2nd and 3rd locations
                Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- TRdf # Add combined TRA to TCR1 related columns 
              }
            }
          } else { # if 1st location is occupied by TRA
            Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- data2[location.i[1],c("TCR1", "cdr3", "cdr3_nt")] # add 1st TRA to TCR1 column
            if (is.na(data2[location.i[2],c("TCR1")])) { # if the 2nd location is occupied by TRB
              if (is.na(data2[location.i[3],c("TCR1")])) { # if the 3rd location is occupied by TRB
                TRdf <- paste(data2[location.i[2],c("TCR2", "cdr3", "cdr3_nt")],data2[location.i[3],c("TCR2", "cdr3", "cdr3_nt")],sep=";") # paste TRB from 2nd and 3rd locations
                Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")] <- TRdf # add TRB combo to TCR2 related  columns
              } else { # if TRA is on 3rd location
                TRdf <- paste(Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")],data2[location.i[3],c("TCR1", "cdr3", "cdr3_nt")],sep=";") # paste TRA from 1st and 3rd locations
                Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- TRdf # add TRA combo to TCR1 related columns 
              }
            } else { # if TRA is on 2nd location
              TRdf <- paste(Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")],data2[location.i[2],c("TCR1", "cdr3", "cdr3_nt")],sep=";") # paste TRA from 1st and 2nd locations
              Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- TRdf # add the combination of 2 TRA chains to the TCR1 slot
              Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")] <- data2[location.i[3],c("TCR2", "cdr3", "cdr3_nt")] # add TRB to the TCR2 related columns 
            }
          }
        } else if (length(location.i) == 1) {
                    chain.i <- data2$chain[location.i]
                    if (chain.i == "TRA"){
                        Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- data2[location.i[1],c("TCR1", "cdr3", "cdr3_nt")]
                    } else {
                        Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")] <- data2[location.i[1],c("TCR2", "cdr3", "cdr3_nt")]
                    }
                }
            }
            Con.df$CTgene <- paste(Con.df$TCR1, Con.df$TCR2, sep="_")
            Con.df$CTnt <- paste(Con.df$cdr3_nt1, Con.df$cdr3_nt2, sep="_")
            Con.df$CTaa <- paste(Con.df$cdr3_aa1, Con.df$cdr3_aa2, sep="_")
            Con.df$CTstrict <- paste(Con.df$TCR1, Con.df$cdr3_nt1, Con.df$TCR2, Con.df$cdr3_nt2, sep="_")
            Con.df[Con.df == "NA_NA"] <- NA #remove the na when gene, aa, or nt is called later
            Con.df[Con.df == "NA_NA_NA_NA"] <- NA #remove the na when nt+gene is called later
            data3 <- merge(data2[,-which(names(data2) %in% c("TCR1","TCR2"))], Con.df, by = "barcode")
            data3 <- data3[, c("barcode", "sample", "ID", "TCR1", "cdr3_aa1", "cdr3_nt1", "TCR2", "cdr3_aa2", "cdr3_nt2", "CTgene", "CTnt", "CTaa", "CTstrict")]
            final[[i]] <- data3
        }
    }
    else if (cells == "B") {
        for (i in seq_along(out)) {

            data2 <- out[[i]]
            data2 <- data2 %>%
                mutate(IGKct = ifelse(chain == "IGK", paste(with(data2, interaction(v_gene,  j_gene, c_gene))), NA)) %>%
                mutate(IGLct = ifelse(chain == "IGL", paste(with(data2, interaction(v_gene,  j_gene, c_gene))), NA)) %>%
                mutate(IGHct = ifelse(chain == "IGH", paste(with(data2, interaction(v_gene,  j_gene, d_gene, c_gene))), NA))
            unique_df <- unique(data2$barcode)
            Con.df <- data.frame(matrix(NA, length(unique_df), 7))
            colnames(Con.df) <- c("barcode","IGH", "cdr3_aa1", "cdr3_nt1", "IGLC", "cdr3_aa2", "cdr3_nt2")
            Con.df$barcode <- unique_df
            y <- NULL
            for (y in 1:length(unique_df)){
                barcode.i <- Con.df$barcode[y]
                location.i <- which(barcode.i == data2$barcode)
                if (length(location.i) == 2){
                    if (is.na(data2[location.i[1],c("IGHct")])) {
                        Con.df[y,c("IGLC", "cdr3_nt2", "cdr3_nt2")] <- data2[location.i[1],c("IGLct", "cdr3", "cdr3_nt")]
                        Con.df[y,c("IGH", "cdr3_nt1", "cdr3_nt1")] <- data2[location.i[2],c("IGHct", "cdr3", "cdr3_nt")]
                    } else {
                        Con.df[y,c("IGH", "cdr3_aa1", "cdr3_nt1")] <- data2[location.i[1],c("IGHct", "cdr3", "cdr3_nt")]
                        Con.df[y,c("IGLC", "cdr3_aa2", "cdr3_nt2")] <- data2[location.i[2],c("IGLct", "cdr3", "cdr3_nt")]
                    }
                } else if (length(location.i) == 1) {
                    chain.i <- data2$chain[location.i]
                    if (chain.i == "IGH"){
                        Con.df[y,c("IGH", "cdr3_aa1", "cdr3_nt1")] <- data2[location.i[1],c("IGHct", "cdr3", "cdr3_nt")]
                    } else {
                        Con.df[y,c("IGLC", "cdr3_aa2", "cdr3_nt2")] <- data2[location.i[2],c("IGLct", "cdr3", "cdr3_nt")]
                    }
                }
            }
            Con.df$CTgene <- paste(Con.df$IGH, Con.df$IGLC, sep="_")
            Con.df$CTnt <- paste(Con.df$cdr3_nt_h, on.df$cdr3_nt_l, sep="_")
            Con.df$CTaa <- paste(on.df$cdr3_h, Con.df$cdr3_l, sep="_")
            Con.df$CTstrict <- paste(Con.df$IGH, Con.df$cdr3_nt_h, Con.df$IGLC, Con.df$cdr3_nt_l, sep="_")
            Con.df[Con.df == "NA_NA"] <- NA #remove the na when gene, aa, or nt is called later
            Con.df[Con.df == "NA_NA_NA_NA"] <- NA #remove the na when nt+gene is called later
            data3 <- merge(data2, Con.df[,-which(names(Con.df) %in% c("IGH","IGLC"))], by = "barcode")
            data3 <- data3[, c("barcode", "sample", "ID", "IGH", "cdr3_aa1", "cdr3_nt1", "IGLC", "cdr3_aa2", "cdr3_nt2", "CTgene", "CTnt", "CTaa", "CTstrict")]
            final[[i]] <- data3
        }
    }
    names <- NULL
    for (i in seq_along(samples)) {
        c <- paste(samples[i], "_", ID[i], sep="")
        names <- c(names, c)}
    names(final) <- names
  
    for (i in seq_along(final)) {
    final[[i]] <- final[[i]][!duplicated(final[[i]]$barcode),]
  }

    return(final)
}
