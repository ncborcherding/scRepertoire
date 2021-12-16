
#' Diversity indices for single-cell RNA-seq
#'
#' @description This function utilizes the Startrac R package derived from 
#' \href{https://pubmed.ncbi.nlm.nih.gov/30479382/}{PMID: 30479382}
#' Required to run the function, the "type" variable needs to include the 
#' difference in where the cells were derived. The output of this function 
#' will produce 3 indices: expa (clonal expansion), migra 
#' (cross-tissue migration), and trans (state transition). In order 
#' to understand the underlying analyses of the outputs please 
#' read and cite the linked manuscript. 
#' 
#' @examples
#' \dontrun{ 
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' #Getting a sample of a Seurat object
#' screp_example <- get(data("screp_example"))
#' screp_example <- combineExpression(combined, screp_example)
#' 
#' #Using occupiedscRepertoire()
#' StartracDiversity(screp_example, type = "Type", sample = "Patient", by = "overall")
#' }
#' 
#' @param sc The seurat or SCE object to visualize after combineExpression(). 
#' For SCE objects, the cluster variable must be in the meta data under 
#' "cluster".
#' @param type The column header in the meta data that gives the where the 
#' cells were derived from, not the patient sample IDs
#' @param sample The column header corresponding to individual samples 
#' or patients. 
#' @param by Method to subset the indices by either overall 
#' (across all samples) or by specific group
#' @param exportTable Returns the data frame used for forming the graph
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
#' @return ggplot object of Startrac diversity metrics
StartracDiversity <- function(sc,
                            type = "Type",
                            sample = NULL,
                            by = "overall", 
                            exportTable = FALSE) {
    meta <- grabMeta(sc)
    colnames(meta)[ncol(meta)] <- "majorCluster"
    meta$clone.status <- ifelse(meta$Frequency > 1, "Yes", "No")
    if (is.null(sample)) {
        stop("Must Add the sample information in order to make 
             the StarTrac calculations")
    } else {
        processed <- data.frame(rownames(meta), meta$CTstrict, meta$clone.status, 
                                meta[,sample], meta[,"majorCluster"], 
                                meta[,type], stringsAsFactors = FALSE)
        colnames(processed) <- c("Cell_Name", "clone.id", "clone.status", "patient", 
                                 "majorCluster", "loc")
    }
    processed <- na.omit(processed)
    indices <- suppressWarnings(Startrac.run(processed, proj = "total", verbose = FALSE))
    indices <- data.frame(indices@cluster.data)
    if (by == "overall") {
        indices <- subset(indices, aid != "total")
        
        melted <- melt(indices)
        values <- as.character(unique(melted$majorCluster))
        values2 <- quiet(dput(values))
        melted$majorCluster <- factor(melted$majorCluster, levels = values2)
        
        plot <- ggplot(melted, aes(x=majorCluster, y=value)) +
            geom_boxplot(aes(fill = majorCluster), outlier.alpha = 0) +
            facet_grid(variable ~.) +
            theme_classic() +
            ylab("Index Score") +
            guides(fill="none") +
            theme(axis.title.x = element_blank())
        
    } else {
        indices <- subset(indices, aid == by)
        plot <- ggplot(melted, aes(x=majorCluster, y=value)) +
            geom_point(aes(fill = by)) +
            facet_grid(variable ~.) +
            theme_classic() +
            ylab("Index Score") +
            guides(fill="none") +
            theme(axis.title.x = element_blank())
    }
    if (exportTable == TRUE) { 
        return(indices) 
        } 
    return(plot)
}




#' The Startrac Class
#'
#' The Startrac object store the data for tcr-based T cell dynamics analyis. The slots contained 
#' in Startrac object are listed below:
#' @slot aid character. aid of the object, used for identification of the object. 
#' For example, patient id. default: "AID"
#' @slot cell.data data.frame. Each line for a cell, and these columns as 
#' required: `Cell_Name`, `clone.id`, `patient`, `majorCluster`, `loc`
#' @slot cell.perm.data object. list of `Startrac`` objects constructed from 
#' permutated cell data
#' @slot clonotype.data data.frame. Each line for a clonotype; contain the 
#' clonotype level indexes information
#' @slot cluster.data data.frame. Each line for a cluster; contain the cluster 
#' level indexes information
#' @slot pIndex.migr data.frame. Each line for a cluster; pairwise migration 
#' index between the two locations indicated in the column name.
#' @slot pIndex.tran data.frame. Each line for a cluster; pairwise transition 
#' index betwwen the two major clusters indicated by the row name and column name.
#' @slot cluster.sig.data data.frame. Each line for a cluster; contains the 
#' p values of cluster indices.
#' @slot pIndex.sig.migr data.frame. Each line for a cluster; contains the 
#' p values of pairwise migration indices.
#' @slot pIndex.sig.tran data.frame. Each line for a cluster; contains the 
#' p values of pairwise transition indices.
#' @slot clonotype.dist.loc matrix. Each line for a clonotype and describe 
#' the cells distribution among the locations.
#' @slot clonotype.dist.cluster matrix. Each line for a clonotype and 
#' describe the cells distribution among the clusters.
#' @slot clust.size array. Number of cells of each major cluster.
#' @slot patient.size array. Number of cells of each patient. 
#' @slot clone.size array. Number of cells of each clone.
#' @slot clone2patient array. Mapping from patient id to clone id.
#' @name Startrac
#' @rdname Startrac
#' @aliases Startrac-class
#' @return method definition for runing startrac
Startrac <- setClass("Startrac",
                    slots = c(aid = "character",
                        cell.data = "data.frame",
                        cell.perm.data = "list",
                        clonotype.data = "data.frame",
                        cluster.data = "data.frame",
                        cluster.sig.data = "data.frame",
                        pIndex.migr = "data.frame",
                        pIndex.tran = "data.frame",
                        pIndex.sig.migr = "data.frame",
                        pIndex.sig.tran = "data.frame",
                        clonotype.dist.loc = "matrix",
                        clonotype.dist.cluster = "matrix",
                        clust.size = "array",
                        patient.size = "array",
                        clone.size = "array",
                        clone2patient = "array"))


setValidity("Startrac",
            function(object) {
                msg <- NULL
                if(!is.data.frame(object@cell.data) || 
                    !all(c("Cell_Name","clone.id","patient","majorCluster","loc") %in% colnames(object@cell.data))){
                    msg <- sprintf("cell.data must be data.frame and contain these columns: Cell_Name, clone.id, 
                                patient, majorCluster, loc")
                }
                if (is.null(msg)) TRUE
                else msg
            }
)


#' show method for Startrac
#
#' @param object A Startrac object
#' @name show
#' @aliases show,Startrac-method 
#' @docType methods
#' @keywords internal
#' @return method for show
setMethod("show",
            signature = "Startrac",
            definition = function(object) {
                cat(sprintf("An object of class %s, aid %s, %d cells from %d clonotypes\n",
                    class(object),
                    object@aid,
                    nrow(object@cell.data),
                    nrow(object@clonotype.data)))
            invisible(x = NULL)
        }
)

#' initialize method for Startrac
#
#' @importFrom plyr llply ldply
#' @param .Object A Startrac object
#' @param cell.data data.frame contains the input data
#' @param aid character analysis id
#' @param n.perm integer number of permutation will be performed. If NULL, no permutation. (default: NULL)
#' @param cores number of core to be used. Passed to doParallel registerDoParallel. default: NULL.
#' @name initialize
#' @aliases initialize,Startrac-method
#' @docType methods
#' @rdname initialize-methods
#' @return an object of class \code{Startrac}
#' @keywords internal
setMethod("initialize",
        signature = "Startrac",
        definition = function(.Object, cell.data, aid="AID",n.perm=NULL,cores=NULL){
            .Object@aid <- aid
            .Object@cell.data <- as.data.frame(cell.data)
            validObject(.Object)
            .Object@clonotype.dist.loc <- unclass(table(.Object@cell.data[,c("clone.id","loc")]))
            .Object@clonotype.dist.cluster <- unclass(table(.Object@cell.data[,c("clone.id","majorCluster")]))
            .Object@clust.size <- unclass(table(.Object@cell.data$majorCluster))
            .Object@patient.size <- unclass(table(.Object@cell.data$patient))
            .Object@clone.size <- unclass(sort(table(.Object@cell.data$clone.id),decreasing = TRUE))
            .clone2patient <- unique(.Object@cell.data[,c("patient","clone.id")])
            .Object@clone2patient <- as.array(structure(.clone2patient$patient, 
                                                names=.clone2patient$clone.id))
            .Object@cell.perm.data <- list()
            if(!is.null(n.perm)){
                registerDoParallel(if(is.null(cores)) (detectCores()-2) else cores)
                .Object@cell.perm.data <- llply(seq_len(n.perm),function(i){
                    perm.cell.data <- ldply(unique(.Object@cell.data$patient),function(pp){
                    .dat <- subset(.Object@cell.data,patient==pp)
                    .dat$clone.id <- .dat$clone.id[sample(nrow(.dat))]
                    return(.dat)
            })
            ###perm.cell.data$clone.id <- perm.cell.data$clone.id[sample(nrow(perm.cell.data))]
            new("Startrac",perm.cell.data,aid=sprintf("perm%06d",i))
            },.progress = "none",.parallel=TRUE)
        }
            return(.Object)
    }
)

#' Calculate cluster level indices
#
#' @name calIndex
#' @aliases calIndex calIndex,Startrac-method
#'
#' @importFrom plyr llply
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @param object A Startrac object
#' @param n.perm integer number of permutation will be performed. If NULL, no permutation. (default: NULL)
#' @param cores number of core to be used. Passed to doParallel registerDoParallel. default: NULL.
#' @param normEntropy logical; whether normalize migration and transition index. default: FALSE.
#' @return an object of class \code{Startrac}
#' @keywords internal
Startrac.calIndex <- function(object,cores,n.perm,normEntropy){
    .entropy <- mcol.entropy(object@clonotype.dist.cluster)
    .entropy.max <- log2(colSums(object@clonotype.dist.cluster > 0))
    object@cluster.data <- data.frame("aid"=object@aid,
                            "majorCluster"=colnames(object@clonotype.dist.cluster),
                            "expa"=1-.entropy/.entropy.max,
                            stringsAsFactors = FALSE)
    if(normEntropy){
        .entropy.migr.max <- log2(ncol(object@clonotype.dist.loc))
        .entropy.tran.max <- log2(ncol(object@clonotype.dist.cluster))
    }else{
        .entropy.migr.max <- 1
        .entropy.tran.max <- 1
    }
    object@clonotype.data <- 
      data.frame("clone.id"=rownames(object@clonotype.dist.loc),
                  "migr"=mrow.entropy(object@clonotype.dist.loc)/.entropy.migr.max,
                  "tran"=mrow.entropy(object@clonotype.dist.cluster)/.entropy.tran.max)
    weights.mtx <- sweep(object@clonotype.dist.cluster,2,colSums(object@clonotype.dist.cluster),"/")
    index.mtx <- t(weights.mtx) %*% (as.matrix(object@clonotype.data[,c("migr","tran")]))
    object@cluster.data <- cbind(object@cluster.data,index.mtx)
    if(!is.null(n.perm)){
        registerDoParallel(if(is.null(cores)) (detectCores()-2) else cores)
        object@cell.perm.data <- llply(object@cell.perm.data,function(x){
            calIndex(x,cores=1,normEntropy=normEntropy)
        },.progress = "none",.parallel=TRUE)
    }
    return(object)
}


setGeneric("calIndex", function(object,cores=NULL,n.perm=NULL,normEntropy=FALSE) standardGeneric("calIndex"))

#' @rdname calIndex
#' @aliases calIndex
#' @keywords internal
setMethod("calIndex", signature = "Startrac", definition = Startrac.calIndex)


#' Calculate pairwise indices
#
#' @name pIndex
#' @aliases pIndex pIndex,Startrac-method
#'
#' @importFrom reshape2 dcast
#' @importFrom plyr ldply adply
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom utils combn
#' @param object A Startrac object
#' @param cores number of core to be used. Passed to doParallel registerDoParallel. default: NULL.
#' @param n.perm integer number of permutation will be performed. If NULL, no permutation. (default: NULL)
#' @return an object of class \code{Startrac}
#' @keywords internal
Startrac.pIndex <- function(object,cores,n.perm){
    clone.dist.loc.majorCluster <- table(object@cell.data[,c("majorCluster","clone.id","loc")])
    withCallingHandlers({
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
    },warning=function(w) {
        if(grepl("... may be used in an incorrect context:",conditionMessage(w)))
            invokeRestart("muffleWarning")
    })
    if(!is.null(cls.migr.index.df) && nrow(cls.migr.index.df)>0){
        cls.migr.index.df$crossLoc <- sprintf("%s-%s",cls.migr.index.df$V1,cls.migr.index.df$V2)
        object@pIndex.migr <- dcast(cls.migr.index.df,majorCluster ~ crossLoc,value.var = "pIndex.migr")
        object@pIndex.migr <- cbind(data.frame(aid=object@aid,stringsAsFactors = FALSE),object@pIndex.migr)
    }else{
        object@pIndex.migr <- data.frame()
    }
    if(!is.null(ncol(object@clonotype.dist.cluster)) && ncol(object@clonotype.dist.cluster)>=2){
        comb.cls <- expand.grid(colnames(object@clonotype.dist.cluster),
                                colnames(object@clonotype.dist.cluster),stringsAsFactors = FALSE)
        comb.cls <- comb.cls[comb.cls[,1]!=comb.cls[,2],]
        withCallingHandlers({
            cls.tran.index.df <- adply(comb.cls,1,function(x,object){
                dat.block <- object@clonotype.dist.cluster[,c(x[[1]],x[[2]])]
                dat.block.clone.index <- mrow.entropy(dat.block)
                dat.block.clone.index[is.na(dat.block.clone.index)] <- 0
                data.frame(pIndex.tran= t(rowSums(dat.block)/sum(dat.block)) %*% dat.block.clone.index)
            },object=object,.progress = "none",.parallel=FALSE)
        },warning=function(w) {
            if(grepl("... may be used in an incorrect context:",conditionMessage(w)))
                invokeRestart("muffleWarning")
        })
        object@pIndex.tran <- dcast(cls.tran.index.df,Var2~Var1,value.var = "pIndex.tran")
        colnames(object@pIndex.tran)[1] <- "majorCluster"
        object@pIndex.tran <- cbind(data.frame(aid=object@aid,stringsAsFactors = FALSE),object@pIndex.tran)
        
        if(!is.null(n.perm)){
            registerDoParallel(if(is.null(cores)) (detectCores()-2)  else cores)
            object@cell.perm.data <- llply(object@cell.perm.data,function(x){
                pIndex(x,n.perm=NULL)
            },.progress = "none",.parallel=TRUE)
        }
    }else{
        object@pIndex.tran <- data.frame()
        if(!is.null(n.perm)){
        }
    }
    return(object)  
}


setGeneric("pIndex", function(object,cores=NULL,n.perm=NULL) standardGeneric("pIndex"))

#' @rdname pIndex
#' @aliases pIndex
#' @keywords internal
setMethod("pIndex", signature = "Startrac", definition = Startrac.pIndex)




#' Get the p value given one Startrac object and a list of Startrac objects from permutation data
#
#' @name getSig
#' 
#' @importFrom plyr laply
#' @importFrom doParallel registerDoParallel
#' @importFrom reshape2 melt
#' @param obj A Startrac object
#' @param obj.perm A list of Startrac objects from permutation data 
#' @keywords internal
#' @return an object of class \code{Startrac}
Startrac.getSig <- function(obj,obj.perm)
{
    .get.empirical.p <- function(obj,obj.perm,slot.name)
    {
        if(nrow(slot(obj,slot.name))==0){
            return(data.frame())
        }
        a.index <- melt(slot(obj,slot.name),id.vars=c("aid","majorCluster"),variable.name="index")
        if(!is.null(obj.perm)){
            perm.mtx <- t(laply(obj.perm,function(x){
                vv <- melt(slot(x,slot.name),id.vars=c("aid","majorCluster"),variable.name="index")
                vv.mtx <- as.matrix(vv[,"value",drop=FALSE])
                rownames(vv.mtx) <- sprintf("%s.%s",vv[["majorCluster"]],vv[["index"]])
                colnames(vv.mtx) <- x@aid
                return(vv.mtx)
            }))
            stopifnot(all(sprintf("%s.%s",a.index$majorCluster,a.index$index)==rownames(perm.mtx)))
            a.index$p.value <- sapply(seq_along(a.index$value),function(i){
                v.rnd <- perm.mtx[i,,drop=FALSE]
                sum(v.rnd >= a.index$value[i])/length(v.rnd)
            })
        }else{
            a.index$p.value <- NA
        }
        return(a.index)
    }
    obj@cluster.sig.data <- .get.empirical.p(obj,obj.perm,"cluster.data")
    obj@pIndex.sig.migr <- .get.empirical.p(obj,obj.perm,"pIndex.migr")
    obj@pIndex.sig.tran <- .get.empirical.p(obj,obj.perm,"pIndex.tran")
    return(obj)
}


setGeneric("getSig", function(obj,obj.perm=NULL) standardGeneric("getSig"))

#' @rdname getSig
#' @aliases getSig
#' @keywords internal
setMethod("getSig", signature = "Startrac", definition = Startrac.getSig)



#' The StartracOUt Class
#'
#' Object store the result of Startrac.run:
#' @slot proj character. identification of the object. For example, patient id. default: "AID"
#' @slot cluster.data data.frame. Each line for a cluster; contain the cluster level indexes information
#' @slot pIndex.migr data.frame. Each line for a cluster; pairwise migration index between the two locations indicated in the column name.
#' @slot pIndex.tran data.frame. Each line for a cluster; pairwise transition index betwwen the two major clusters indicated by the row name and column name.
#' @slot cluster.sig.data data.frame. Each line for a cluster; contains the p values of cluster indices.
#' @slot pIndex.sig.migr data.frame. Each line for a cluster; contains the p values of pairwise migration indices.
#' @slot pIndex.sig.tran data.frame. Each line for a cluster; contains the p values of pairwise transition indices.
#' @slot objects list. other objects
#' @name StartracOut
#' @rdname StartracOut
#' @aliases StartracOut-class
#' @exportClass StartracOut
#' @keywords internal
StartracOut <- setClass("StartracOut",
                        slots = c(proj = "character",
                        cluster.data = "data.frame",
                        cluster.sig.data = "data.frame",
                        pIndex.migr = "data.frame",
                        pIndex.tran = "data.frame",
                        pIndex.sig.migr = "data.frame",
                        pIndex.sig.tran = "data.frame",
                        objects = "list"))

#' initialize method for StartracOut
#
#' @param .Object A StartracOut object
#' @param proj character analysis id
#' @aliases initialize,StartracOut-method
#' @docType methods
#' @keywords internal
#' @return an object of class \code{StartracOut}
setMethod("initialize",
            signature = "StartracOut",
            definition = function(.Object, proj="AID"){
                .Object@proj <- proj
            return(.Object)
            }
)


#' show method for StartracOut
#
#' @importFrom utils head
#' @param object A StartracOut object
#' @aliases show,StartracOut-method
#' @docType methods
#' @keywords internal
setMethod("show",
            signature = "StartracOut",
            definition = function(object) {
            cat(sprintf("An object of class %s, proj %s:\n",
                    class(object),
                    object@proj))
            cat("head of the clusters' index:\n")
            print(head(object@cluster.data))
            cat("head of the pairwise migration index:\n")
            print(head(object@pIndex.migr))
            cat("head of the pairwise transition index:\n")
            print(head(object@pIndex.tran[,seq_len(min(5,ncol(object@pIndex.tran)))]))
            invisible(x = NULL)
        }
)



#' dispaly message with time stamp
#' @param msg characters; message to display
#' @keywords internal
#' @return Estimate of sys.time
loginfo <- function(msg) {
    timestamp <- sprintf("%s", Sys.time())
    msg <- paste0("[",timestamp, "] ", msg,"\n")
    cat(msg)
}

#' entropy of each row of the input matrix
#' @param x matrix;
#' @keywords internal
#' @return row entropy
mrow.entropy <- function(x)
{
    freqs <- sweep(x,1,rowSums(x),"/")
    H = - rowSums(ifelse(freqs>0,freqs* log2(freqs),0))
    return(H)
}

#' entropy of each column of the input matrix
#' @param x matrix;
#' @keywords internal
#' @return  column entropy
mcol.entropy <- function(x)
{
    freqs <- sweep(x,2,colSums(x),"/")
    H = - colSums(ifelse(freqs>0,freqs* log2(freqs),0))
    return(H)
}

#' warpper function for Startrac analysis
#' @importFrom reshape2 dcast
#' @importFrom plyr ldply adply llply
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom methods new slot slot<- new validObject
#' @param cell.data data.frame. Each line for a cell, and these columns as required: `Cell_Name`, `clone.id`, `patient`, `majorCluster`, `loc`
#' @param proj character. String used to annotate the project.
#' @param cores integer. number of core to be used. default: NULL.
#' @param n.perm integer. number of permutation will be performed. If NULL, no permutation. (default: NULL)
#' @param verbose logical. wheter return intermediate result (some Startrac objects) 
#' @details run the Startrac pipeline
#' @keywords internal
#' @return an list contains data.frame elements "cluster.data","pIndex.migr" and "pIndex.tran"

Startrac.run <- function(cell.data, proj="CRC", cores=NULL,n.perm=NULL,verbose=FALSE)
{
    ##tic("obj.proj")
    loginfo("initialize Startrac ...")
    obj.proj <- new("Startrac",cell.data,aid=proj,n.perm=n.perm,cores=cores)
    loginfo("calculate startrac index ...")
    obj.proj <- calIndex(obj.proj,cores=cores,n.perm=n.perm)
    loginfo("calculate pairwise index ...")
    obj.proj <- pIndex(obj.proj,cores=cores,n.perm=n.perm)
    if(!is.null(n.perm)){ 
        loginfo("get the significance")
        obj.proj <- getSig(obj.proj,obj.proj@cell.perm.data) 
    }else{
        obj.proj <- getSig(obj.proj,NULL) 
    }
    ##toc()
    
    obj.list <- NULL
    if(length(obj.proj@patient.size)>1)
    {
        loginfo("calculate indices of each patient ...")
        patient.vec <- names(obj.proj@patient.size[obj.proj@patient.size > 30])
        #cl <- makeCluster(if(is.null(cores)) (detectCores()-2) else cores)
        #registerDoParallel(cl)
        withCallingHandlers({
            obj.list <- llply(patient.vec,function(pid,cell.data){
                obj <- new("Startrac",subset(cell.data,patient==pid),aid=pid)
                obj <- calIndex(obj)
                obj <- pIndex(obj,cores=1)
                obj <- getSig(obj,NULL)
                obj
            },cell.data=cell.data,.progress = "none",.parallel=FALSE)
        },warning=function(w) {
            if(grepl("... may be used in an incorrect context:",conditionMessage(w)))
                ### strange bug, see https://github.com/hadley/plyr/issues/203
                invokeRestart("muffleWarning")
        })
        #stopCluster(cl)
    }
    loginfo("collect result")
    ret <- new("StartracOut",proj=proj)
    ## cluster index
    ret.slot.names <- c("cluster.data","pIndex.migr","pIndex.tran",
                        "cluster.sig.data","pIndex.sig.migr","pIndex.sig.tran")
    #  if(!is.null(n.perm)){ 
    #    ret.slot.names <- c(ret.slot.names,
    #                        c("cluster.sig.data","pIndex.sig.migr","pIndex.sig.tran")) 
    #  }
    for(v in ret.slot.names)
    {
        slot(ret, v) <- slot(obj.proj,v)
        if(!is.null(obj.list)){
            slot(ret, v) <- rbind(slot(ret, v),ldply(obj.list,function(obj){
                slot(obj,v)
            }))
        }
    }
    if(verbose){
        ret@objects <- c(obj.proj,obj.list)
    }
    loginfo("return")
    return(ret)
}
