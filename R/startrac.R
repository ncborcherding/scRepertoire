StartracDiversity <- function(seurat,
                              type = "Type",
                              sample = NULL,
                              by = c("overall")) {
    meta <- data.frame(seurat[[]], Idents(seurat), stringsAsFactors = F)
    colnames(meta)[ncol(meta)] <- "majorCluster"
    meta$clone.status <- ifelse(meta$Frequency > 1, "Yes", "No")
    if (is.null(sample)) {
        stop("Must Add the sample information in order to make the StarTrac calculations")
    } else {
        processed <- data.frame(rownames(meta), meta$CTstrict, meta$clone.status, meta[,sample], meta[,"majorCluster"], meta[,type], stringsAsFactors = F)
        colnames(processed) <- c("Cell_Name", "clone.id", "clone.status", "patient", "majorCluster", "loc")
    }
    processed <- na.omit(processed)
    indices <- suppressWarnings(Startrac.run(processed, proj = "total", verbose = F))
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
            guides(fill=F) +
            theme(axis.title.x = element_blank())
        
    } else {
        indices <- subset(indices, aid == by)
        plot <- ggplot(melted, aes(x=majorCluster, y=value)) +
            geom_point(aes(fill = by)) +
            facet_grid(variable ~.) +
            theme_classic() +
            ylab("Index Score") +
            guides(fill=F) +
            theme(axis.title.x = element_blank())
    }
    
    return(plot)
    
}




#' The Startrac Class
#'
#' The Startrac object store the data for tcr-based T cell dynamics analyis. The slots contained in Startrac object are listed below:
#' @slot aid character. aid of the object, used for identification of the object. For example, patient id. default: "AID"
#' @slot cell.data data.frame. Each line for a cell, and these columns as required: `Cell_Name`, `clone.id`, `patient`, `majorCluster`, `loc`
#' @slot cell.perm.data object. list of `Startrac`` objects constructed from permutated cell data
#' @slot clonotype.data data.frame. Each line for a clonotype; contain the clonotype level indexes information
#' @slot cluster.data data.frame. Each line for a cluster; contain the cluster level indexes information
#' @slot pIndex.migr data.frame. Each line for a cluster; pairwise migration index between the two locations indicated in the column name.
#' @slot pIndex.tran data.frame. Each line for a cluster; pairwise transition index betwwen the two major clusters indicated by the row name and column name.
#' @slot cluster.sig.data data.frame. Each line for a cluster; contains the p values of cluster indices.
#' @slot pIndex.sig.migr data.frame. Each line for a cluster; contains the p values of pairwise migration indices.
#' @slot pIndex.sig.tran data.frame. Each line for a cluster; contains the p values of pairwise transition indices.
#' @slot clonotype.dist.loc matrix. Each line for a clonotype and describe the cells distribution among the locations.
#' @slot clonotype.dist.cluster matrix. Each line for a clonotype and describe the cells distribution among the clusters.
#' @slot clust.size array. Number of cells of each major cluster.
#' @slot patient.size array. Number of cells of each patient. 
#' @slot clone.size array. Number of cells of each clone.
#' @slot clone2patient array. Mapping from patient id to clone id.
#' @name Startrac
#' @rdname Startrac
#' @aliases Startrac-class
#' @exportClass Startrac
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
                    msg <- sprintf("cell.data must be data.frame and contain these columns: Cell_Name, clone.id, patient, majorCluster, loc")
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
#' @param cores number of core to be used. Passed to doParallel::registerDoParallel. default: NULL.
#' @name initialize
#' @aliases initialize,Startrac-method
#' @docType methods
#' @rdname initialize-methods
#' @return an object of class \code{Startrac}
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
              .Object@clone.size <- unclass(sort(table(.Object@cell.data$clone.id),decreasing = T))
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
                  },.progress = "none",.parallel=T)
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
#' @param cores number of core to be used. Passed to doParallel::registerDoParallel. default: NULL.
#' @param normEntropy logical; whether normalize migration and transition index. default: FALSE.
#' @return an object of class \code{Startrac}
Startrac.calIndex <- function(object,cores,n.perm,normEntropy)
{
    ### cluster level expansion index (STARTRAC-expa)
    #### Todo: special case: number of clonotype is 1, i.e. sum(x>0)==1
    .entropy <- mcol.entropy(object@clonotype.dist.cluster)
    .entropy.max <- log2(colSums(object@clonotype.dist.cluster > 0))
    object@cluster.data <- data.frame("aid"=object@aid,
                                      "majorCluster"=colnames(object@clonotype.dist.cluster),
                                      "expa"=1-.entropy/.entropy.max,
                                      stringsAsFactors = F)
    ### clone level migration and transition index
    if(normEntropy){
        .entropy.migr.max <- log2(ncol(object@clonotype.dist.loc))
        .entropy.tran.max <- log2(ncol(object@clonotype.dist.cluster))
    }else{
        .entropy.migr.max <- 1
        .entropy.tran.max <- 1
    }
    object@clonotype.data <- data.frame("clone.id"=rownames(object@clonotype.dist.loc),
                                        "migr"=mrow.entropy(object@clonotype.dist.loc)/.entropy.migr.max,
                                        "tran"=mrow.entropy(object@clonotype.dist.cluster)/.entropy.tran.max)
    ### cluster level migration index (STARTRAC-migr) and transition index (STARTRAC-tran)
    weights.mtx <- sweep(object@clonotype.dist.cluster,2,colSums(object@clonotype.dist.cluster),"/")
    index.mtx <- t(weights.mtx) %*% (as.matrix(object@clonotype.data[,c("migr","tran")]))
    object@cluster.data <- cbind(object@cluster.data,index.mtx)
    if(!is.null(n.perm)){
        #cl <- makeCluster(if(is.null(cores)) (detectCores()-2) else cores)
        #registerDoParallel(cl)
        registerDoParallel(if(is.null(cores)) (detectCores()-2) else cores)
        object@cell.perm.data <- llply(object@cell.perm.data,function(x){
            calIndex(x,cores=1,normEntropy=normEntropy)
        },.progress = "none",.parallel=T)
        #stopCluster(cl)
        
    }
    return(object)
}

#' @export
setGeneric("calIndex", function(object,cores=NULL,n.perm=NULL,normEntropy=FALSE) standardGeneric("calIndex"))

#' @rdname calIndex
#' @aliases calIndex
setMethod("calIndex", signature = "Startrac", definition = Startrac.calIndex)


#' Calculate pairwise indices
#
#' @name pIndex
#' @aliases pIndex pIndex,Startrac-method
#'
#' @importFrom data.table dcast
#' @importFrom plyr ldply adply
#' @importFrom utils combn
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @param object A Startrac object
#' @param cores number of core to be used. Passed to doParallel::registerDoParallel. default: NULL.
#' @param n.perm integer number of permutation will be performed. If NULL, no permutation. (default: NULL)
#' @return an object of class \code{Startrac}
Startrac.pIndex <- function(object,cores,n.perm)
{
    ####### index given two cluster or loc
    ## migr 
    clone.dist.loc.majorCluster <- table(object@cell.data[,c("majorCluster","clone.id","loc")])
    
    #tic("migr")
    #clone.dist.loc.majorCluster.debug <<- clone.dist.loc.majorCluster
    withCallingHandlers({
        cls.migr.index.df <- ldply(seq_len(dim(clone.dist.loc.majorCluster)[1]),function(i,clone.dist.loc.majorCluster){
            dat.cls <- clone.dist.loc.majorCluster[i,,]
            i.name <- dimnames(clone.dist.loc.majorCluster)[["majorCluster"]][i]
            dat.cls.index <- NULL
            if(!is.null(ncol(dat.cls)) && ncol(dat.cls)>=2){
                comb.loc <- as.data.frame(t(combn(colnames(dat.cls),2)),stringsAsFactors=F)
                dat.cls.pIndex.migr <- apply(comb.loc,1,function(x){
                    dat.block <- dat.cls[,x]
                    dat.block.clone.index <- mrow.entropy(dat.block)
                    dat.block.clone.index[is.na(dat.block.clone.index)] <- 0
                    t(rowSums(dat.block)/sum(dat.block)) %*% dat.block.clone.index
                })
                dat.cls.index <- cbind(data.frame(majorCluster=rep(i.name,nrow(comb.loc)),
                                                  pIndex.migr=dat.cls.pIndex.migr,
                                                  stringsAsFactors = F),
                                       comb.loc)
            }
            return(dat.cls.index)
        },clone.dist.loc.majorCluster=clone.dist.loc.majorCluster,.progress = "none",.parallel=F)
    },warning=function(w) {
        if(grepl("... may be used in an incorrect context:",conditionMessage(w)))
            ### strange bug, see https://github.com/hadley/plyr/issues/203
            invokeRestart("muffleWarning")
    })
    #toc()
    if(!is.null(cls.migr.index.df) && nrow(cls.migr.index.df)>0){
        cls.migr.index.df$crossLoc <- sprintf("%s-%s",cls.migr.index.df$V1,cls.migr.index.df$V2)
        object@pIndex.migr <- dcast(cls.migr.index.df,majorCluster ~ crossLoc,value.var = "pIndex.migr")
        object@pIndex.migr <- cbind(data.frame(aid=object@aid,stringsAsFactors = F),object@pIndex.migr)
    }else{
        object@pIndex.migr <- data.frame()
    }
    
    ## tran
    if(!is.null(ncol(object@clonotype.dist.cluster)) && ncol(object@clonotype.dist.cluster)>=2){
        ##comb.cls <- as.data.frame(t(combn(colnames(object@clonotype.dist.cluster),2)),stringsAsFactors=F)
        comb.cls <- expand.grid(colnames(object@clonotype.dist.cluster),
                                colnames(object@clonotype.dist.cluster),stringsAsFactors = F)
        comb.cls <- comb.cls[comb.cls[,1]!=comb.cls[,2],]
        
        #tic("tran")
        withCallingHandlers({
            cls.tran.index.df <- adply(comb.cls,1,function(x,object){
                dat.block <- object@clonotype.dist.cluster[,c(x[[1]],x[[2]])]
                dat.block.clone.index <- mrow.entropy(dat.block)
                dat.block.clone.index[is.na(dat.block.clone.index)] <- 0
                data.frame(pIndex.tran= t(rowSums(dat.block)/sum(dat.block)) %*% dat.block.clone.index)
            },object=object,.progress = "none",.parallel=F)
        },warning=function(w) {
            if(grepl("... may be used in an incorrect context:",conditionMessage(w)))
                ### strange bug, see https://github.com/hadley/plyr/issues/203
                invokeRestart("muffleWarning")
        })
        #toc()
        
        object@pIndex.tran <- dcast(cls.tran.index.df,Var2~Var1,value.var = "pIndex.tran")
        colnames(object@pIndex.tran)[1] <- "majorCluster"
        object@pIndex.tran <- cbind(data.frame(aid=object@aid,stringsAsFactors = F),object@pIndex.tran)
        
        if(!is.null(n.perm)){
            #cl <- makeCluster(if(is.null(cores)) (detectCores()-2)  else cores)
            #registerDoParallel(cl)
            registerDoParallel(if(is.null(cores)) (detectCores()-2)  else cores)
            object@cell.perm.data <- llply(object@cell.perm.data,function(x){
                pIndex(x,n.perm=NULL)
            },.progress = "none",.parallel=T)
            #stopCluster(cl)
        }
    }else{
        object@pIndex.tran <- data.frame()
        if(!is.null(n.perm)){
            ### nothing needed to do
        }
    }
    return(object)  
}

#' @export
setGeneric("pIndex", function(object,cores=NULL,n.perm=NULL) standardGeneric("pIndex"))

#' @rdname pIndex
#' @aliases pIndex
setMethod("pIndex", signature = "Startrac", definition = Startrac.pIndex)




#' Get the p value given one Startrac object and a list of Startrac objects from permutation data
#
#' @name getSig
#' @aliases getSig getSig,Startrac-method
#' 
#' @importFrom plyr laply
#' @importFrom doParallel registerDoParallel
#' @importFrom data.table melt
#' @param obj A Startrac object
#' @param obj.perm A list of Startrac objects from permutation data 
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
                vv.mtx <- as.matrix(vv[,"value",drop=F])
                rownames(vv.mtx) <- sprintf("%s.%s",vv[["majorCluster"]],vv[["index"]])
                colnames(vv.mtx) <- x@aid
                return(vv.mtx)
            }))
            stopifnot(all(sprintf("%s.%s",a.index$majorCluster,a.index$index)==rownames(perm.mtx)))
            a.index$p.value <- sapply(seq_along(a.index$value),function(i){
                v.rnd <- perm.mtx[i,,drop=F]
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

#' @export
setGeneric("getSig", function(obj,obj.perm=NULL) standardGeneric("getSig"))

#' @rdname getSig
#' @aliases getSig
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
              print(head(object@pIndex.tran[,1:min(5,ncol(object@pIndex.tran))]))
              invisible(x = NULL)
          }
)

#' plot the indexes
#
#' @name plot
#' @aliases plot plot,StartracOut-method
#' 
#' @importFrom plyr laply
#' @importFrom doParallel registerDoParallel
#' @importFrom data.table melt as.data.table melt
#' @importFrom ggpubr ggbarplot ggboxplot
#' @importFrom ggplot2 facet_wrap theme element_text aes geom_text
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize colorRamp2
#' @importFrom grDevices colorRampPalette
#' @param obj A object of StartracOut
#' @param index.type one of "cluster.all", "pairwise.migr", "pairwise.tran". (default:"cluster.all")
#' @param byPatient logical. plot indexes of each patient (default: FALSE)
#' @return a ggplot2 object or Heatmap-class object
StartracOut.plot <- function(obj,index.type,byPatient)
{
    if(index.type=="cluster.all"){
        if(byPatient){
            p <- ggboxplot(as.data.table(obj@cluster.sig.data)[aid!=obj@proj,][order(majorCluster),],
                           x="majorCluster",y="value",palette = "npg",
                           color = "index", add = "point", outlier.colour=NULL) +
                facet_wrap(~index,ncol=1,scales = "free_y") +
                theme(axis.text.x=element_text(angle = 60,hjust = 1))
        }else{
            dat.plot <- as.data.table(obj@cluster.sig.data)[aid==obj@proj,]
            dat.plot$p.value.label <- ""
            dat.plot$p.value.label[dat.plot$p.value < 0.05] <- "*"
            dat.plot$p.value.label[dat.plot$p.value < 0.01] <- "**"
            dat.plot$p.value.label[dat.plot$p.value < 0.001] <- "***"
            p <- ggbarplot(dat.plot[order(majorCluster),],
                           x="majorCluster",y="value",palette = "npg",fill = "index") +
                facet_wrap(~index,ncol=1,scales = "free_y") +
                coord_cartesian(clip="off") +
                theme(axis.text.x=element_text(angle = 60,hjust = 1),strip.background = element_blank())
            if(!all(is.na(dat.plot$p.value))){
                p <- p + geom_text(aes(label=p.value.label,y=value),size=5)
            }
        }
        
    }else if(index.type=="pairwise.migr"){
        if(byPatient){
            p <- ggboxplot(as.data.table(obj@pIndex.sig.migr)[aid!=obj@proj,][order(majorCluster),],
                           x="majorCluster",y="value",palette = "npg",
                           color = "index", add = "point", outlier.colour=NULL) +
                facet_wrap(~index,ncol=1,scales = "free_y") +
                theme(axis.text.x=element_text(angle = 60,hjust = 1))      
        }else{
            dat.plot <- as.data.table(obj@pIndex.sig.migr)[aid==obj@proj,]
            dat.plot$p.value.label <- ""
            dat.plot$p.value.label[dat.plot$p.value < 0.05] <- "*"
            dat.plot$p.value.label[dat.plot$p.value < 0.01] <- "**"
            dat.plot$p.value.label[dat.plot$p.value < 0.001] <- "***"
            p <- ggbarplot(dat.plot[order(majorCluster),],
                           x="majorCluster",y="value",palette = "npg",fill = "index") +
                facet_wrap(~index,ncol=1,scales = "free_y") +
                coord_cartesian(clip="off") +
                theme(axis.text.x=element_text(angle = 60,hjust = 1),strip.background = element_blank())
            if(!all(is.na(dat.plot$p.value))){
                p <- p + geom_text(aes(label=p.value.label,y=value),size=5)
            }
        }
    }else if(index.type=="pairwise.tran"){
        dat.plot <- as.matrix(subset(obj@pIndex.tran,aid==obj@proj)[,c(-1,-2)])
        rownames(dat.plot) <- subset(obj@pIndex.tran,aid==obj@proj)[,2]
        dat.plot[is.na(dat.plot)] <- 0
        yrange <- pretty(dat.plot)
        col.heat <- colorRamp2(seq(0,max(yrange),length=15),
                               colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(15),
                               space = "LAB")
        p <- Heatmap(dat.plot,name="pIndex.tran",col = col.heat)
    }
    return(p)
}

#' @export
setGeneric("plot", function(obj,index.type="cluster.all",byPatient=F) standardGeneric("plot"))

#' @rdname plot
#' @aliases plot
setMethod("plot", signature = "StartracOut", definition = StartracOut.plot)



#' dispaly message with time stamp
#' @param msg characters; message to display
#' @export
loginfo <- function(msg) {
    timestamp <- sprintf("%s", Sys.time())
    msg <- paste0("[",timestamp, "] ", msg,"\n")
    cat(msg)
}

#' entropy of each row of the input matrix
#' @param x matrix;
mrow.entropy <- function(x)
{
    freqs <- sweep(x,1,rowSums(x),"/")
    H = - rowSums(ifelse(freqs>0,freqs* log2(freqs),0))
    return(H)
}

#' entropy of each column of the input matrix
#' @param x matrix;
mcol.entropy <- function(x)
{
    freqs <- sweep(x,2,colSums(x),"/")
    H = - colSums(ifelse(freqs>0,freqs* log2(freqs),0))
    return(H)
}

#' warpper function for Startrac analysis
#' @importFrom data.table dcast
#' @importFrom plyr ldply adply llply
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom methods new slot
#' @importFrom methods slot<-
#' @param cell.data data.frame. Each line for a cell, and these columns as required: `Cell_Name`, `clone.id`, `patient`, `majorCluster`, `loc`
#' @param proj character. String used to annotate the project.
#' @param cores integer. number of core to be used. default: NULL.
#' @param n.perm integer. number of permutation will be performed. If NULL, no permutation. (default: NULL)
#' @param verbose logical. wheter return intermediate result (some Startrac objects) 
#' @details run the Startrac pipeline
#' @return an list contains data.frame elements "cluster.data","pIndex.migr" and "pIndex.tran"
#' @export
#' @examples 
#' library("Startrac")
#' dat.file <- system.file("extdata/example.cloneDat.Zhang2018.txt",package = "Startrac")
#' in.dat <- read.table(dat.file,stringsAsFactors = FALSE,head=TRUE)
#' out <- Startrac.run(in.dat, proj="CRC", cores=2,verbose=FALSE)
#' 
Startrac.run <- function(cell.data, proj="CRC", cores=NULL,n.perm=NULL,verbose=F)
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
                require("Startrac")
                obj <- new("Startrac",subset(cell.data,patient==pid),aid=pid)
                obj <- calIndex(obj)
                obj <- pIndex(obj,cores=1)
                obj <- getSig(obj,NULL)
                obj
            },cell.data=cell.data,.progress = "none",.parallel=F)
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

#' calculate Startrac.dist (tissue distribution preference)
#' @import data.table
#' @importFrom plyr aaply
#' @importFrom stats chisq.test
#' @param dat.tb data.frame. Each line for a cell, and these columns as required: `majorCluster`, `loc`
#' @param byPatient logical. whether calculate the index for each patient. (default: FALSE)
#' @param colname.cluster character. which column specify the cluster (default: "majorCluster")
#' @param colname.patient character. which column specify the patient  (default: "patient")
#' @param colname.tissue character. which column specify the tissue  (default: "loc")
#' @details calculate Startrac.dist (tissue distribution preference) which is based on Chisquare test.
#' @return an array full of R_{o/e}
#' @export
calTissueDist <- function(dat.tb,byPatient=F,colname.cluster="majorCluster",
                          colname.patient="patient",colname.tissue="loc")
{
    if(byPatient==F){
        N.o <- table(dat.tb[[colname.cluster]],dat.tb[[colname.tissue]])
        res.chisq <- chisq.test(N.o)
        R.oe <- (res.chisq$observed)/(res.chisq$expected)
    }else{
        N.o.byPatient <- table(dat.tb[[colname.patient]],
                               dat.tb[[cluster.colname]], dat.tb[[colname.tissue]])
        R.oe <- aaply(N.o.byPatient,1,function(x){
            res.chisq <- chisq.test(x)
            return((res.chisq$observed)/(res.chisq$expected))
        })
    }
    return(R.oe)
}