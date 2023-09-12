#' Startrac-based diversity indices for single-cell RNA-seq
#'
#' @description This function utilizes the Startrac approach derived from 
#' \href{https://pubmed.ncbi.nlm.nih.gov/30479382/}{PMID: 30479382}
#' Required to run the function, the "type" variable needs to include the 
#' difference in where the cells were derived. The output of this function 
#' will produce 3 indices: expa (clonal expansion), migra 
#' (cross-tissue migration), and trans (state transition). In order 
#' to understand the underlying analyses of the outputs please 
#' read and cite the linked manuscript. 
#' 
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' #Getting a sample of a Seurat object
#' scRep_example  <- get(data("scRep_example "))
#' scRep_example  <- combineExpression(combined, scRep_example )
#' scRep_example$Patient <- substring(scRep_example$orig.ident,1,3)
#' scRep_example$Type <- substring(scRep_example$orig.ident,4,4) 
#' 
#' #Using occupiedscRepertoire()
#' StartracDiversity(scRep_example , type = "Type", sample = "Patient", by = "overall")
#'
#' 
#' @param sc The seurat or SCE object to visualize after combineExpression(). 
#' For SCE objects, the cluster variable must be in the meta data under 
#' "cluster".
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param type The variable in the meta data that provides tissue type
#' @param group.by The variable in the meta data to group by, often samples
#' @param metric Calculate the "standard" values across the group.by variable
#'  or "pairwise" values
#' @param exportTable Returns the data frame used for forming the graph
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
#' @return ggplot object of Startrac diversity metrics
StartracDiversity <- function(sc,
                              cloneCall = "strict", 
                              chain = "both",
                              type = "Type",
                              group.by = NULL, 
                              metric = "standard",
                              exportTable = FALSE) {
    
    cloneCall <- theCall(cloneCall)
    df <- grabMeta(sc)
    barcodes <- rownames(df)
    colnames(df)[ncol(df)] <- "majorCluster"
    
    if (is.null(group.by)) {
       if("orig.ident" %!in% colnames(df)) {
         stop("Please select a group.by variable")
       }
       group.by <- "orig.ident"
    }
    group.levels = unique(dat[,group.by])
    
    cloneCall <- theCall(cloneCall)
    if (chain != "both") {
      df <- off.the.chain(df, chain, cloneCall)
    }

    df <- df %>%
              group_by(df[,group.by], df[,cloneCall]) %>%
              dplyr::mutate(n = n()) %>%
              as.data.frame()
                
    rownames(df) <- barcodes
    remove.pos <- which(df[,cloneCall] %in% c("", NA))
    if (length(remove.pos) > 0) {
      df <- df[-remove.pos,]
    }
    df[,"clone.status"] <- ifelse(df[,"n"] > 1, "Yes", "No")
    
    if (is.null(sample)) {
        stop("Must Add the sample information in order to make 
             the StarTrac calculations")
    } else {
        processed <- data.frame(rownames(df), df[,cloneCall], df$clone.status, 
                                df[,group.by], df[,"majorCluster"], 
                                df[,type], stringsAsFactors = FALSE)
        colnames(processed) <- c("Cell_Name", "clone.id", "clone.status", "patient", 
                                 "majorCluster", "loc")
    }
    processed[processed == "NA"] <- NA
    processed <- na.omit(processed)
    indexFunc <- switch(metric,
                        "standard" = .calIndex,
                        "pairwise"  = .pIndex,
                        stop("Invalid metric provided"))
    
    mat.list <- list()
    for(i in seq_along(group.levels)) {
        mat.list[[i]] <- as.data.frame(indexFunc(processed[processed[,4] == group.levels[i],]))
    }
    names(mat.list) <- group.levels
    if(metric == "standard") {
      mat <- bind_rows(mat.list, .id = "group")
    } else if (metric == "pairwise") {
      tran.list <- bind_rows(mat)
    }
    
    if (exportTable) { 
      return(mat) 
    } 
    
        
    mat_melt <- melt(mat)
    values <- as.character(unique(mat_melt$majorCluster))
    values2 <- quiet(dput(values))
    mat_melt$majorCluster <- factor(mat_melt$majorCluster, levels = values2)
        
    plot <- ggplot(mat_melt, aes(x=majorCluster, y=value)) +
            geom_boxplot(aes(fill = majorCluster), outlier.alpha = 0) +
            facet_grid(variable ~.) +
            theme_classic() +
            ylab("Index Score") +
            guides(fill="none") +
            theme(axis.title.x = element_blank())
        
    
   
    return(plot)
}


#' Calculate cluster level indices
#' @importFrom plyr llply
#' @keywords internal
.calIndex <- function(processed){
    clonotype.dist.cluster <- table(processed[,c("clone.id","majorCluster")])
    clonotype.dist.loc <- unclass(table(processed[,c("clone.id","loc")]))
    .entropy <- mcol.entropy(clonotype.dist.cluster)
    .entropy.max <- log2(colSums(clonotype.dist.cluster > 0))
    expa <- 1-.entropy/.entropy.max

    .entropy.migr.max <- 1
    .entropy.tran.max <- 1
  
    
    clonotype.data <- 
      data.frame("clone.id"=rownames(clonotype.dist.loc),
                  "migr"=mrow.entropy(clonotype.dist.loc)/.entropy.migr.max,
                  "tran"=mrow.entropy(clonotype.dist.cluster)/.entropy.tran.max)
    weights.mtx <- sweep(clonotype.dist.cluster,2,colSums(clonotype.dist.cluster),"/")
    calIndex.matrix <- t(weights.mtx) %*% (as.matrix(clonotype.data[,c("migr","tran")]))
    calIndex.matrix <- cbind(calIndex.matrix, expa = expa)
    calIndex.matrix[is.nan(calIndex.matrix)] <- NA
    return(calIndex.matrix)
}



#' Calculate pairwise indices
#' @importFrom reshape2 dcast
#' @importFrom plyr ldply adply
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom utils combn
.pIndex <- function(processed){
    clone.dist.loc.majorCluster <- table(processed[,c("majorCluster","clone.id","loc")])

        cls.migr.index.df <- ldply(seq_len(dim(clone.dist.loc.majorCluster)[1]),
                                function(i,clone.dist.loc.majorCluster){
            dat.cls <- clone.dist.loc.majorCluster[i,,]
            i.name <- dimnames(clone.dist.loc.majorCluster)[["majorCluster"]][i]
            dat.cls.index <- NULL
            if(!is.null(ncol(dat.cls)) && ncol(dat.cls)>=2){
                comb.loc <- as.data.frame(t(combn(colnames(dat.cls),2)),stringsAsFactors=FALSE)
                dat.cls.pIndex.migr <- apply(comb.loc,1,function(x){
                    dat.block <- dat.cls[,x]
                    dat.block.clone.index <- mrow.entropy(dat.block)
                    dat.block.clone.index[is.na(dat.block.clone.index)] <- 0
                    t(rowSums(dat.block)/sum(dat.block)) %*% dat.block.clone.index
                })
                dat.cls.index <- cbind(data.frame(majorCluster=rep(i.name,nrow(comb.loc)),
                                    pIndex.migr=dat.cls.pIndex.migr,
                                    stringsAsFactors = FALSE),
                                    comb.loc)
            }
            return(dat.cls.index)
        },clone.dist.loc.majorCluster=clone.dist.loc.majorCluster,.progress = "none",.parallel=FALSE)

    if(!is.null(cls.migr.index.df) && nrow(cls.migr.index.df)>0){
        cls.migr.index.df$crossLoc <- sprintf("%s-%s",cls.migr.index.df$V1,cls.migr.index.df$V2)
        pIndex.migr <- dcast(cls.migr.index.df,majorCluster ~ crossLoc,value.var = "pIndex.migr")
    }else{
        pIndex.migr <- data.frame()
    }
    clonotype.dist.cluster <- table(processed[,c("clone.id","majorCluster")])
    if(!is.null(ncol(clonotype.dist.cluster)) && ncol(clonotype.dist.cluster)>=2){
        comb.cls <- expand.grid(colnames(clonotype.dist.cluster),
                                colnames(clonotype.dist.cluster),stringsAsFactors = FALSE)
        comb.cls <- comb.cls[comb.cls[,1]!=comb.cls[,2],]
        cls.tran.index.df <- adply(comb.cls,1,function(x,object){
                dat.block <- clonotype.dist.cluster[,c(x[[1]],x[[2]])]
                dat.block.clone.index <- mrow.entropy(dat.block)
                dat.block.clone.index[is.na(dat.block.clone.index)] <- 0
                data.frame(pIndex.tran= t(rowSums(dat.block)/sum(dat.block)) %*% dat.block.clone.index)
            },object=object,.progress = "none",.parallel=FALSE)
        
        pIndex.tran <- as.data.frame(dcast(cls.tran.index.df,Var2~Var1,value.var = "pIndex.tran"))
        colnames(pIndex.tran)[1] <- "majorCluster"
        
    } else{
        pIndex.tran <- data.frame()
    }
    
    pIndex.matrix <- list("tran" = pIndex.tran, 
                          "migr" = pIndex.migr)
    return(pIndex.matrix)  
}

#' entropy of each row of the input matrix
#' @param x matrix;
#' @keywords internal
#' @return row entropy
mrow.entropy <- function(x)
{
    freqs <- sweep(x,1,rowSums(x),"/")
    H <- - rowSums(ifelse(freqs>0,freqs* log2(freqs),0))
    return(H)
}

#' entropy of each column of the input matrix
#' @param x matrix;
#' @keywords internal
#' @return  column entropy
mcol.entropy <- function(x)
{
    freqs <- sweep(x,2,colSums(x),"/")
    H <- - colSums(ifelse(freqs>0,freqs* log2(freqs),0))
    return(H)
}
