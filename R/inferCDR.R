#' Add portions on the CDR loop based on Vgene usage
#'
#' This function will use the Vgene to add the amino acid sequence 
#' of the CDR loop for given clones. For example, recovering the 
#' CDR1 and CDR2 sequences, will produce a string with 
#' CDR1-CDR2-CDR3 (ex: "SGH.......RS-YFS....ETQ-CASSLTDRTYEQYF"). 
#' Sequences are derived from IMGT-based nomenclature. If a 
#' vgene is not present or does not match IMGT, unrecoverable
#' regions will be denoted with "!" (ex: "!-!-CASSLTDRTYEQYF").
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' TRB_CDR123 <- inferCDR(combined[[1]], 
#'                        chain = "TRB", 
#'                        regions = c("CDR1", "CDR2"), 
#'                        species = "human")
#'
#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param regions CDR loop regions - "FR1", "CDR1", "FR2", "CDR2", and/or "FR3" 
#' @param species The species of the experiment - "human" or "mouse"
#' @importFrom stringr str_split
#' @export
#' @return vector of sequences corresponding with the additional regions

inferCDR <- function(df, 
                     chain = "TRB", 
                     regions = c("CDR1", "CDR2"), 
                     species = "human") {
  
  .sequence.positions[["FR3.F"]] <- c(97:103)
  sco <- is_seurat_object(df) | is_se_object(df)
  df <- .data.wrangle(df, group.by, "CTgene", chain)
  for (x in seq_along(df)) {
    df[[x]] <- .off.the.chain(df[[x]], chain, "CTaa")
    df[[x]] <- .off.the.chain(df[[x]], chain, "CTgene")
  }
  #Get the sequence position for the regions selected
  sequence.pos <- .sequence.positions[grep(paste0(regions, collapse = "|"), names(.sequence.positions))]
  
  dat <- bind_rows(df)
  cdr3.sequences <- dat[,"CTaa"]
  Vgene <- str_split(dat[,"CTgene"], "[.]", simplify = TRUE)[,1]
  
  #Laoding reference
  protein.reference <- .load.fasta(chain, species)
  
  #Matching the tenx Vgene nomenclature with IMGT
  Vgene <- .tenX.V(Vgene, protein.reference)
  
  for (i in seq_len(length(regions)-1)) {
    sequence.pos[[i]] <- c(sequence.pos[[i]], "-")
  }
  
  sequence.pos <- unname(unlist(sequence.pos))
  
  #Loop through the sequences
  sequence.store <- c()
  for(i in seq_len(length(cdr3.sequences))) {
    if (Vgene[i] %in% names(protein.reference)) {
      aa_reference <- protein.reference[[Vgene[i]]]
      aa_sequence <- suppressWarnings(aa_reference[as.numeric(sequence.pos)])
      aa_sequence[is.na(aa_sequence)] <- "-"
      aa_sequence <- paste0(aa_sequence, collapse = "")
      tmp.sequence <- paste0(aa_sequence, "-", cdr3.sequences[i])
    } else {
      aa_sequence <- rep("!", length(regions))
      aa_sequence <- paste0(aa_sequence, collapse = "-")
      tmp.sequence <- paste0(aa_sequence, "-", cdr3.sequences[i])
    }
    sequence.store <- c(sequence.store, tmp.sequence)
  }
  sequence.store[sequence.store == "!-!-NA"] <- NA
  return(sequence.store)
}

#Positions of specific regions
.sequence.positions <- list(FR1.A = c(1:15), 
                            FR1.B = c(16:26),
                            CDR1 = c(27:38), 
                            FR2.C = c(39:46), 
                            FR2.Cp = c(47:55),
                            CDR2 = c(56:65), 
                            FR3.C = c(66:74), 
                            FR3.D = c(75:84),
                            FR3.E = c(85:96), 
                            FR3.F = c(97:104))

#' @importFrom stringr str_remove_all
.tenX.V <- function(Vgene, reference) {
  Vgene <- str_remove_all(Vgene, "NA_")
  #Persistent Adaptive Conversion
  V.dict <- as.data.frame(table(Vgene))
  unmatched.V <- as.vector(V.dict$Vgene[V.dict$Vgene %!in% names(reference)])
  rematched.V <- unmatched.V
  for(i in seq_along(unmatched.V)) {
    tmp.V <- unmatched.V[i]
    tmp.V <- str_sort(names(reference), numeric = TRUE)[grep(tmp.V, str_sort(names(reference), numeric = TRUE))][1]
    if(is.na(tmp.V)) {
      tmp.V <- unmatched.V[i]
      tmp.V <- str_split(tmp.V, "-", simplify = TRUE)[,1]
      tmp.V <- str_sort(names(reference), numeric = TRUE)[grep(tmp.V, str_sort(names(reference), numeric = TRUE))][1]
    }
    rematched.V[i] <- tmp.V
  }
  
  #Replacing Sticky Genes
  names(rematched.V) <- unmatched.V
  matching_indices <- match(Vgene, names(rematched.V))
  Vgene[!is.na(matching_indices)] <- rematched.V[matching_indices[!is.na(matching_indices)]]
  return(Vgene)
}

#Pulling the reference sequences form IMGT
#' @importFrom seqinr read.fasta translate
#' @importFrom stringr str_split str_replace_all
.load.fasta <- function(chain, species) {
  species <- tolower(species)
  species.link <- switch(species,
                          "human" = "Homo_sapiens",
                          "mouse" = "Mus_musculus",
                          stop("Please select either human or mouse for species"))
  
  chain.link <- switch(chain,
                         "TRB" = "TR/TRBV",
                         "TRA" = "TR/TRAV",
                         "TRG" = "TR/TRGV",
                         "TRD" = "TR/TRDV",
                         "Heavy" = "IG/IGHV", 
                         "Light" = c("IG/IGLV", "IG/IGKV"),
                         stop("Please select chain: TRB, TRA, TRG, TRD, Heavy, Light"))
  
  dir <- paste0("https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/",
                species.link, "/", chain.link, ".fasta")
  reference <- lapply(dir, read.fasta)
  if(chain == "Light") {
    reference <- c(reference[[1]], reference[[2]])
  } else {
    reference <- reference[[1]]
  }
  names(reference) <- str_split(names(reference), "[|]", simplify = TRUE)[,2]
  protein.reference <- lapply(reference, function(x) {
    tmp <- seqinr::translate(x)
    tmp <- str_replace_all(tmp, "X", ".")
    tmp
  })
  
}
