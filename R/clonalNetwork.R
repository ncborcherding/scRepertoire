#' Visualize clonal network along reduced dimensions
#'
#' This function generates a network based on clonal 
#' proportions of an indicated identity and then superimposes
#' the network onto a single-cell object dimensional reduction
#' plot. 
#' 
#' @examples

#' \dontrun{
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' #Getting a sample of a Seurat object
#' screp_example <- get(data("screp_example"))
#' 
#' #Using combineExpresion()
#' screp_example <- combineExpression(combined, screp_example)
#' 
#' #Using clonalNetwork()
#' clonalNetwork(screp_example, reduction = "umap",
#'               identity = "cluster")
#' }
#'               
#' @param sc The Seurat or SingleCellExperiment (SCE) after combineExpression().
#' @param reduction The name of the dimensional reduction of the single-cell object
#' @param identity A variable in the meta data to use for the nodes.
#' @param filter.clones Use to select the top n clones (filter.clones = 2000) or 
#' n of clones based on the minimum number of all the comparators (filter.clone = "min").
#' @param filter.identity Display the network for a specific level of the indicated identity
#' @param filter.proportion Remove clonotypes from the network below a specific proportion
#' @param filter.graph Remove the reciprocal edges from the half of the graph, 
#' allowing for cleaner visualization
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param exportTable Exports a table of the data into the global 
#' environment in addition to the visualization
#' @import ggplot2
#' @importFrom stringr str_sort
#' @importFrom igraph graph_from_data_frame V `V<-`
#' @importFrom dplyr %>% group_by select summarize_all count
#' @importFrom tidygraph as_tbl_graph activate
#' @importFrom ggraph ggraph geom_edge_bend  geom_node_point scale_edge_colour_gradientn circle guide_edge_colourbar
#' @export
#' @return ggplot object
#' 
#' 
clonalNetwork <- function(sc, 
                          reduction = "umap",
                          identity = "ident",
                          filter.clones = NULL,
                          filter.identity = NULL,
                          filter.proportion = NULL,
                          filter.graph = FALSE,
                          cloneCall = "strict", 
                          chain = "both", 
                          exportTable = FALSE) {
    
    cloneCall <- theCall(cloneCall)
    meta <- grabMeta(sc)  
    coord <- data.frame(get.coord(sc, reduction), identity = meta[,identity])
    min <- c()
    if (!is.null(filter.clones))  {
      if(filter.clones == "min") {
        meta <- grabMeta(sc)
        id.meta <- split(meta, meta[,identity])
        for (x in seq_along(id.meta)) {
            min.tmp <- length(which(!is.na(unique(id.meta[[x]][,cloneCall]))))
            min <- c(min.tmp, min)
        }
        min <- min(min)
        table <- meta %>%
            group_by(meta[,identity], meta[, cloneCall]) %>%
            count() %>%
            na.omit() %>%
            arrange(desc(n)) 
        runSum <-as.vector(cumsum(table[,3])) 
        table$cumSum <- as.vector(runSum$n)
        cut <- which.min(unlist(abs(table[,4] - min)))
        clones.to.filter <- table[seq_len(cut),1]
    } else if (is.numeric(filter.clones)) {
        table <- meta %>%
          count(meta[, cloneCall]) %>%
          na.omit() %>%
          arrange(desc(n)) 
        runSum <-as.vector(cumsum(table[,2])) 
        table$cumSum <- as.vector(runSum)
        cut <- which.min(unlist(abs(table[,3] - filter.clones)))
        clones.to.filter <- table[seq_len(cut),1]
    }  
    } 
    meta <- grabMeta(sc)    
    if (!is.null(filter.clones)) {
      meta <- meta[meta[,cloneCall] %in% clones.to.filter,]
    }
    clones.duplicated <- na.omit(unique(meta[which(duplicated(meta[,cloneCall])),cloneCall]))
    id <- as.vector(meta[,identity])
    id.names <- id
    names(id.names) <- row.names(meta)
    unique.id <- unique(id)
    
    id.positions <- data.frame(coord)
    colnames(id.positions)[1:2] <- c("x", "y")
    suppressMessages(centers <- id.positions%>% 
        group_by(identity) %>% 
        select(x, y) %>% 
        summarize_all(mean) %>%
        data.frame())
    rownames(centers) <- centers[,1]
    centers <- centers[,-1]
    
    clone.number <- meta[,c(cloneCall, identity)] %>%
        group_by(meta[,identity]) %>% 
        na.omit() %>%
        unique() %>%
        summarise(n = n())
    names <- clone.number[,1]
    clone.number <- unlist(clone.number[,2])
    names(clone.number) <- unlist(names[,1])
    
    total.number <- meta[,c(cloneCall, identity)] %>%
        group_by(meta[,identity]) %>% 
        na.omit() %>%
      summarise(n = n())
    names <- total.number[,1]
    total.number <- unlist(total.number[,2])
    names(total.number) <- unlist(names[,1])
    
    edge.list <- NULL
    for (i in seq_along(clones.duplicated)) {
       pos <- which(meta[,cloneCall] == clones.duplicated[i])
       num <- table(meta[pos,identity])
       num <- num[num > 0]
       if(length(num) == 1) {
         next()
       }
       grid <- expand.grid(names(num),names(num))
       grid <- grid[grid[,1] != grid[,2],]
       for (x in seq_len(nrow(grid))) {
          summary <- c(to = as.vector(grid[x,1]), from = as.vector(grid[x,2]), weight = num[grid[x,2]]/total.number[as.vector(grid[x,1])])
          edge.list <- rbind(edge.list, summary)
       }
    }
    if(is.null(edge.list)) {
      stop("No shared clonotypes between the indicated identity")
    }
    
    edge.list <- data.frame(edge.list)
    colnames(edge.list)[3] <-"weight"
   
    
    if (!is.null(filter.identity)) { 
        col1 <- which(edge.list[,1] == filter.identity)
        col2 <- which(edge.list[,2] == filter.identity)
        edge.list <- edge.list[c(col1,col2),]
    } else if (filter.graph) {
        unique.id <- str_sort(unique.id, numeric = TRUE)
        edge.list <- edge.list[edge.list[,1] %in% unique.id[seq_len(length(unique.id)/2)],]
    }
    if(!is.null(filter.proportion)) {
        edge.list <-edge.list[edge.list[,3] > filter.proportion,]
    }
    edge.list1 <- edge.list %>%
        group_by(to, from) %>%
        summarise(weight = sum(as.numeric(weight)))
    graph <- graph_from_data_frame(edge.list1)
    clone.number <- clone.number[match(names(V(graph)), names(clone.number))]
    V(graph)$size <- unname(clone.number)
    centers <- centers[rownames(centers) %in% names(V(graph)),]

    plot1 <- ggraph(graph, layout = centers[match(names(V(graph)), rownames(centers)),]) + 
        geom_point(data = coord, aes(x = coord[,1], y = coord[,2], color = identity)) + 
        geom_edge_bend(aes(edge_color = as.numeric(weight)), alpha = 0.7, width = 1,
                       arrow = arrow(length = unit(4, 'mm')), 
                      end_cap = circle(3, 'mm'), 
                      angle_calc = "across", 
                      check_overlap = TRUE, 
                      strength = 0.75) + 
        geom_node_point(aes(size = size)) + 
        ylab(paste0(reduction, "_2")) + 
        xlab(paste0(reduction, "_1")) + 
        guides(color = "none") + 
        scale_edge_colour_gradientn(colors  = rev(colorblind_vector(13)), trans = "log10") + 
        labs(size = "Unique Clones", edge_color = "Relative Proportion of \nClones of Starting Node") + 
        theme(
            panel.background = element_blank(),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"), 
            legend.key=element_blank()
        )
    if (exportTable == TRUE) { return(edge.list1)}
    return(plot1)
}

