#' Hierarchical clustering of clonotypes on clonotype size and 
#' Jensen-Shannon divergence
#'
#' This function produces a hierarchical clustering of clonotypes by sample 
#' using the Jensen-Shannon distance and discrete gamma-GPD spliced threshold 
#' model in 
#' \href{https://bioconductor.org/packages/devel/bioc/html/powerTCR.html}{powerTCR R package}.
#' Please read and cite PMID: 30485278 if using the function for analyses. 
#' 
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' clonalSizeDistribution(combined, cloneCall = "strict", method="ward.D2")
#'
#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param cloneCall How to call the clonotype - VDJC gene (gene),
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param threshold Numerical vector containing the thresholds 
#' the grid search was performed over.
#' @param method The clustering parameter for the dendrogram.
#' @param group.by The column header used for grouping.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any hcl.pals()
#' @import ggplot2 
#' @importFrom dplyr bind_rows
#' @importFrom ggdendro dendro_data segment label
#' @importFrom stats hclust optim pgamma
#' @importFrom stats as.dist
#' @export
#' @return ggplot dendrogram of the clone size distribution
#' @author Hillary Koch

clonalSizeDistribution <- function(df,
                                   cloneCall ="strict", 
                                   chain = "both", 
                                   method = "ward.D2", 
                                   threshold = 1, 
                                   group.by = NULL,
                                   split.by = NULL, 
                                   exportTable = FALSE, 
                                   palette = "inferno") {
  x <- xend <- yend <- mpg_div_hp <- NULL
  cloneCall <- .theCall(cloneCall)
  df <- .data.wrangle(df, split.by, cloneCall, chain)
  
  if(!is.null(group.by)) {
    df <- .groupList(df, group.by)
  }
  data <- bind_rows(df)
  unique_df <- unique(data[,cloneCall])
  
  #Forming data frame to store values
  Con.df <- data.frame(matrix(NA, length(unique_df), length(df)))
  Con.df <- data.frame(unique_df, Con.df, stringsAsFactors = FALSE)
  colnames(Con.df)[1] <- "clonotype"
  for (i in seq_along(df)) {
    data <- df[[i]]
    data <- data.frame(table(data[,cloneCall]), 
                       stringsAsFactors = FALSE)
    colnames(data) <- c(cloneCall, "Freq")
    for (y in seq_along(unique_df)){
      clonotype.y <- Con.df$clonotype[y]
      location.y <- which(clonotype.y == data[,cloneCall])
      Con.df[y,i+1] <- data[location.y[1],"Freq"]
    }
  }
  colnames(Con.df)[2:(length(df)+1)] <- names(df)
  Con.df[is.na(Con.df)] <- 0
  list <- list()
  for (i in seq_along(df)) {
    list[[i]] <- suppressWarnings(.fdiscgammagpd(Con.df[,i+1], useq = threshold))
  }
  names(list) <- names(df)
  grid <- 0:10000
  distances <- .get_distances(list, grid, modelType="Spliced")
  mat_melt <- dendro_data(hclust(as.dist(distances), method = method), type = "rectangle")
  
  #Plotting
  plot <- ggplot() + 
            geom_segment(data = segment(mat_melt), 
                         aes(x = x, y = y, xend = xend, yend = yend)) +
            geom_text(data = label(mat_melt), 
                              aes(x = x, y = -0.02, label = label, hjust = 0), size = 4) +
            geom_point(data = label(mat_melt), 
                              aes(x = x, y = -0.01, color = as.factor(label)), size = 2) + 
            coord_flip() +
            scale_y_reverse(expand = c(0.2, 0)) + 
            scale_color_manual(values = .colorizer(palette, nrow(label(mat_melt)))) + 
            theme_classic() + 
            guides(color = "none") + 
            theme(axis.title = element_blank(), 
                  axis.ticks.y = element_blank(), 
                  axis.text.y = element_blank()) 
  
  if (exportTable) { 
    return(distances) 
  }
  return(plot)
}

#####################################################
#Functions to run the fdiscgammagpd and get_distances
#####################################################
#This is to reduce dependency footprint of scRepertoire
#and ensure compatibility going forward.
#' @importFrom evmix fgpd
.fdiscgammagpd <- function(x, useq, shift = NULL, pvector=NULL,
                          std.err = TRUE, method = "Nelder-Mead", ...){
  if(!is(x, "numeric")){
    stop("x must be numeric.")
  }
  
  if(!is.null(shift)){
    if(!is(shift, "numeric")){
      stop("shift must be numeric.")
    }
    if(shift != round(shift)){
      stop("shift must be an integer.")
    }
  }
  
  if(!is(useq, "numeric")){
    stop("useq must be numeric.")
  }
  
  if(any(x != round(x))){
    stop("all elements in x must be integers.")
  }
  
  if(any(useq != round(useq))){
    stop("all elements in useq must be integers.")
  }
  
  if(!is.null(pvector) & !(length(pvector) == 5)){
    stop("pvector must contain 5 elements.")
  }
  
  if(!(is.logical(std.err))){
    stop("std.err must be TRUE or FALSE.")
  }
  
  if(is.null(shift)){
    shift <- min(x)
  }
  
  if (is.null(pvector)) {
    pvector <- rep(NA,5)
    s <- log(mean(x+0.5))-mean(log(x+0.5))
    k <- (3-s + sqrt((s-3)^2 + 24*s))/12/s
    pvector[1] <- k
    pvector[2] <- k/mean(x)
    pvector[3] <- as.vector(quantile(x, 0.9))
    
    xu <- x[x>=pvector[3]]
    initfgpd <- fgpd(xu, min(xu)-10^(-5))
    pvector[4] <- initfgpd$mle[1]
    pvector[5] <- initfgpd$mle[2]
  }
  
  bulk <- lapply(seq_along(useq),
                 function(idx,x,useq) x < useq[idx], x=x, useq=useq)
  tail <- lapply(seq_along(useq),
                 function(idx,x,useq) x >= useq[idx], x=x, useq=useq)
  phiu <- lapply(seq_along(useq),
                 function(idx,tail) mean(tail[[idx]]), tail=tail)
  
  gammfit <- list()
  gpdfit <- list()
  nllhu <- rep(NA, length(useq))
  for(i in seq_along(useq)){
    gammfit[[i]] <- tryCatch(expr = .fdiscgamma(pvector[1:2],x[bulk[[i]]],
                                               useq[i],
                                               phiu[[i]],
                                               shift,
                                               method="Nelder-Mead"),
                             error = function(err) NA)
    gpdfit[[i]] <- tryCatch(expr = .fdiscgpd(pvector[4:5],
                                            x[tail[[i]]],
                                            useq[i],
                                            phiu[[i]],
                                            method="Nelder-Mead"),
                            error = function(err) {
                              pvec3 <- as.vector(quantile(x,1-phiu[[i]]))
                              xu <- x[x>=pvec3]
                              initfgpd.adj <-
                                fgpd(x, min(xu)-10^(-5))
                              pvec4 <- initfgpd.adj$mle[1]
                              pvec5 <- initfgpd.adj$mle[2]
                              tryCatch(expr = .fdiscgpd(c(pvec4,pvec5),
                                                       x[tail[[i]]],
                                                       useq[i],
                                                       phiu[[i]],
                                                       method="Nelder-Mead"),
                                       error = function(err2) NA)
                            })
    nllhu[i] <- tryCatch(expr = gammfit[[i]]$value + gpdfit[[i]]$value,
                         error = function(err) NA)
  }
  
  bestfit <- which.min(nllhu)
  fit.out <- list(gammfit[[bestfit]], gpdfit[[bestfit]])
  names(fit.out) <- c("bulk","tail")
  mle <- c(mean(x >= useq[bestfit]),
           exp(fit.out$bulk$par),
           useq[bestfit],
           exp(fit.out$tail$par[1]),
           fit.out$tail$par[2])
  names(mle) <- c("phi","shape","rate","thresh","sigma","xi")
  if(std.err){
    H <- fit.out$bulk$hessian %>% rbind(matrix(rep(0,4),nrow = 2)) %>%
      cbind(rbind(matrix(rep(0,4),nrow = 2),fit.out$tail$hessian))
    fisherInf <- tryCatch(expr = solve(H), error = function(err) NA)
    out <- list(x = as.vector(x), shift = shift, init = as.vector(pvector),
                useq = useq, nllhuseq = nllhu,
                optim = fit.out, nllh = nllhu[bestfit],
                mle=mle, fisherInformation = fisherInf)
  } else{
    out <- list(x = as.vector(x), shift = shift, init = as.vector(pvector),
                useq = useq, nllhuseq = nllhu,
                optim = fit.out, nllh = nllhu[bestfit],
                mle=mle)
  }
  out
}
##-----------------------------------------------------------------------------
## Density, distribution, and quantile functions, random number generation
## for discrete truncated gamma and discrete gpd
##-----------------------------------------------------------------------------

.ddiscgamma <- function(x, shape, rate, thresh, phiu, shift = 0, log = FALSE){
  if(any(x != floor(x))){
    stop("x must be an integer")
  }
  
  out <- rep(0, length(x))
  
  up <- pgamma(x+1-shift, shape=shape, rate=rate)
  down <- pgamma(x-shift, shape=shape, rate=rate)
  
  if(!log){
    b <- pgamma(thresh-shift, shape=shape, rate=rate)
    out[x < thresh] <- ((1-phiu)*(up-down)/b)[x < thresh]
  } else{
    b <- pgamma(thresh-shift, shape=shape, rate=rate, log.p = TRUE)
    out[x < thresh] <- (log(1-phiu)+log(up-down) - b)[x < thresh]
  }
  out
}

.pdiscgamma <- function(q, shape, rate, thresh, phiu, shift = 0){
  probs <- .ddiscgamma(0:q, shape, rate, thresh, phiu, shift)
  sum(probs)
}

#' @importFrom truncdist qtrunc
.qdiscgamma <- function(p, shape, rate, thresh, phiu, shift = 0){
  qtrunc(p/(1-phiu), spec = "gamma", a=0,
         b=thresh-shift, shape=shape,rate=rate) %>% floor + shift
}

#' @importFrom truncdist rtrunc
.rdiscgamma <- function(n, shape, rate, thresh, shift = 0){
  rtrunc(n, spec = "gamma", a=0,
         b=thresh-shift, shape=shape, rate=rate) %>% floor + shift
}

#' @importFrom evmix pgpd
.ddiscgpd <- function(x, thresh, sigma, xi, phiu, log = FALSE){
  up <- pgpd(x+1, u=thresh, sigmau=sigma, xi=xi)
  down <- pgpd(x, u=thresh, sigmau=sigma, xi=xi)
  
  if(!log){
    phiu*(up-down)
  } else{
    log(phiu) + log(up-down)
  }
}

pdiscgpd <- function(q, thresh, sigma, xi, phiu){
  probs <- .ddiscgpd(thresh:q, thresh, sigma, xi, phiu)
  sum(probs)
}

#' @importFrom evmix qgpd
.qdiscgpd <- function(p, thresh, sigma, xi, phiu){
  qgpd(p/phiu, u=thresh, sigmau = sigma, xi=xi) %>% floor
}

#' @importFrom evmix rgpd
.rdiscgpd <- function(n, thresh, sigma, xi){
  rgpd(n, u=thresh, sigmau=sigma, xi=xi) %>% floor
}


##----------------------------------------------------------------------------------
## negative log likelihood and parameter estimation functions
## for discrete truncated gamma and discrete gpd
##----------------------------------------------------------------------------------

.discgammanll <- function(param, dat, thresh, phiu, shift=0){
  shape <- exp(param[1])
  rate <- exp(param[2])
  
  if(any(dat > thresh-1)){ warning("data must be less than the threshold") }
  
  ll <- log(.ddiscgamma(dat, shape, rate, thresh, phiu, shift))
  
  sum(-ll)
}

.discgpdnll <- function(param, dat, thresh, phiu){
  sigma <- exp(param[1])
  xi <- param[2]
  ll <- log(.ddiscgpd(dat, thresh, sigma, xi, phiu))
  
  sum(-ll)
}

.fdiscgamma <- function(param, dat, thresh, phiu, shift = 0, method, ...){
  opt <- optim(log(param), .discgammanll, dat=dat, thresh=thresh,
               phiu=phiu, shift=shift, method=method, hessian = TRUE, ...)
  opt
}

.fdiscgpd <- function(param, dat, thresh, phiu, method, ...){
  opt <- optim(c(log(param[1]),param[2]), .discgpdnll, dat=dat, thresh=thresh,
               phiu=phiu, method=method, hessian = TRUE, ...)
  opt
}

.get_distances <- function(fits, grid, modelType = "Spliced"){
  if(!is(grid, "numeric")){
    stop("grid must be numeric.")
  }
  if(any(grid != round(grid))){
    stop("all elements in grid must be integers.")
  }
  if(!(modelType %in% c("Spliced", "Desponds"))){
    stop("modelType must be either \"Spliced\" or \"Desponds\".")
  }
  
  distances <- matrix(rep(0, length(fits)^2), nrow = length(fits))
  if(!is.null(names(fits))){
    rownames(distances) <- colnames(distances) <- names(fits)
  }
  
  for(i in seq_len((length(fits)-1))){
    for(j in (i+1):length(fits)){
      distances[i,j] <- .JS_dist(fits[[i]],
                                fits[[j]],
                                grid,
                                modelType = modelType)
    }
  }
  distances <- distances + t(distances)
  distances
}

.JS_dist <- function(fit1, fit2, grid, modelType = "Spliced"){
  if(!(modelType %in% c("Spliced", "Desponds"))){
    stop("modelType must be either \"Spliced\" or \"Desponds\".")
  }
  
  if(modelType == "Spliced"){
    if(!all(c("x", "init", "useq", "nllhuseq", "nllh",
              "optim", "mle") %in% names(fit1))){
      stop("\"fit1\" is not of the correct structure. It must be a model
                 fit from fdiscgammagpd.")
    }
    if(!all(c("x", "init", "useq", "nllhuseq", "nllh",
              "optim", "mle") %in% names(fit2))){
      stop("\"fit2\" is not of the correct structure. It must be a model
                 fit from fdiscgammagpd.")
    }
    shiftp <- fit1$shift
    shiftq <- fit2$shift
    phip <- fit1$mle['phi']
    phiq <- fit2$mle['phi']
    shapep <- fit1$mle['shape']
    shapeq <- fit2$mle['shape']
    ratep <- fit1$mle['rate']
    rateq <- fit2$mle['rate']
    threshp <- fit1$mle['thresh']
    threshq <- fit2$mle['thresh']
    sigmap <- fit1$mle['sigma']
    sigmaq <- fit2$mle['sigma']
    xip <- fit1$mle['xi']
    xiq <- fit2$mle['xi']
    
    out <- .JS_spliced(grid, shiftp, shiftq, phip, phiq, shapep, shapeq,
                      ratep, rateq, threshp, threshq, sigmap, sigmaq,
                      xip, xiq)
  } else if(modelType == "Desponds"){
    if(!all(c("min.KS", "Cmin", "powerlaw.exponent",
              "pareto.alpha") == names(fit1))){
      stop("\"fit1\" is not of the correct structure. It must be a model
                 fit from fdesponds.")
    }
    if(!all(c("min.KS", "Cmin", "powerlaw.exponent",
              "pareto.alpha") == names(fit2))){
      stop("\"fit2\" is not of the correct structure. It must be a model
                 fit from fdesponds.")
    }
    Cminp <- fit1['Cmin']
    Cminq <- fit2['Cmin']
    alphap <- fit1['pareto.alpha']
    alphaq <- fit2['pareto.alpha']
    out <- .JS_desponds(grid, Cminp, Cminq, alphap, alphaq)
  }
  out
}

#' @importFrom evmix dgpd
.JS_spliced <- function(grid, shiftp, shiftq, phip, phiq, shapep, shapeq, ratep,
                       rateq, threshp, threshq, sigmap, sigmaq, xip, xiq){
  if(!is(grid, "numeric")){
    stop("grid must be numeric.")
  }
  
  if(any(grid != round(grid))){
    stop("all elements in grid must be integers.")
  }
  
  if(any(!is(c(shiftp, shiftq, phip, phiq, shapep, shapeq,
               ratep, rateq, threshp, threshq,
               sigmap, sigmaq, xip, xiq), "numeric"))){
    stop("shiftp, shiftq, phip, phiq, shapep, shapeq, ratep, rateq,
              threshp, threshq, sigmap, sigmaq, xip, and xiq must be numeric.")
  }
  
  if(shiftp != round(shiftp) | shiftq != round(shiftq)){
    stop("shiftp and shiftq must be integers.")
  }
  
  if(any(c(shapep, shapeq, ratep, rateq, sigmap, sigmaq) <= 0)){
    stop("shapep, shapeq, ratep, rateq, sigmap, and sigmaq must be 
             greater than 0.")
  }
  
  if(any(c(phip, phiq) > 1) | any(c(phip, phiq) < 0)){
    stop("phip and phiq must be in [0,1].")
  }
  
  if(ratep <= 0 | rateq <= 0){
    stop("ratep and rateq must be greater than 0.")
  }
  
  if(shapep <= 0 | shapeq <= 0){
    stop("shapep and shapeq must be greater than 0.")
  }
  
  if(threshp != round(threshp) | threshq != round(threshq)){
    stop("threshp and threshq must be integers.")
  }
  
  K <- max(grid)
  
  P <- .ddiscgammagpd(min(grid):K, shape = shapep, rate = ratep,
                     u=threshp, sigma = sigmap,
                     xi = xip, phiu = phip, shift=shiftp,
                     log = FALSE)
  adjp <- which(P == 0)
  if(length(adjp) != 0){
    P[adjp] <- dgpd(adjp+0.5, u=threshp,
                           sigmau = sigmap, xi = xip, phiu = phip)
  }
  
  Q <- .ddiscgammagpd(min(grid):K, shape = shapeq, rate = rateq,
                     u=threshq, sigma = sigmaq,
                     xi = xiq, phiu = phiq, shift=shiftq,
                     log = FALSE)
  adjq <- which(Q == 0)
  if(length(adjq) != 0){
    Q[adjq] <- dgpd(adjq+0.5, u=threshq,
                           sigmau = sigmaq, xi = xiq, phiu = phiq)
  }
  
  M <- 0.5*(P+Q)
  pzero <- which(P == 0)
  qzero <- which(Q == 0)
  
  sum1 <- sum2 <- rep(NA, length(grid))
  sum1 <- P*(log(P) - log(M))
  sum2 <- Q*(log(Q) - log(M))
  
  if(length(intersect(pzero, qzero)) != 0){
    sum1[intersect(pzero, qzero)] <- 0
    sum2[intersect(pzero, qzero)] <- 0
  }
  if(length(setdiff(pzero, qzero)) != 0){
    sum1[setdiff(pzero, qzero)] <- 0
  }
  if(length(setdiff(qzero, pzero)) != 0){
    sum2[setdiff(qzero, pzero)] <- 0
  }
  
  out <- sqrt(0.5*(sum(sum1) + sum(sum2)))
  out
}


#' @importFrom cubature adaptIntegrate
.JS_desponds <- function(grid, Cminp, Cminq, alphap, alphaq){
  if(!is(grid, "numeric")){
    stop("grid must be numeric.")
  }
  
  if(any(grid != round(grid))){
    stop("all elements in grid must be integers.")
  }
  
  if(any(!is(c(Cminp, Cminq, alphap, alphaq), "numeric"))){
    stop("Cminp, Cminq, alphap, and alphaq must be numeric.")
  }
  
  if(Cminp != round(Cminp) | Cminq != round(Cminq)){
    stop("Cminp and Cminq must be integers.")
  }
  
  if(alphap <= 0 | alphaq <= 0){
    stop("alphap and alphaq must be greater than 0.")
  }
  
  lower = min(grid)
  upper = max(grid)
  
  out <- adaptIntegrate(.eval_desponds,
                        lowerLimit = lower, upperLimit = upper,
                        Cminp = Cminp, Cminq = Cminq,
                        alphap = alphap, alphaq = alphaq)$integral %>% sqrt
  out
}

.ddiscgammagpd <- function(x, fit=NULL, shape, rate, u, sigma, xi,
                          phiu=NULL, shift = 0, log = FALSE){
  if(!is.null(fit)){
    if(!all(c("x", "shift", "init", "useq", "nllhuseq", "nllh",
              "optim", "mle") %in% names(fit))){
      stop("\"fit\" is not of the correct structure. It must be one model
                 fit from fdiscgammagpd.")
    }
    phiu <- fit$mle['phi']
    shape <- fit$mle['shape']
    rate <- fit$mle['rate']
    u <- fit$mle['thresh']
    sigma <- fit$mle['sigma']
    xi <- fit$mle['xi']
    shift <- fit$shift
  }
  if(!is(x, "numeric")){
    stop("x must be numeric.")
  }
  
  if(!is(shift, "numeric")){
    stop("shift must be numeric.")
  }
  
  if(any(x != floor(x))){
    stop("x must be an integer")
  }
  
  if(shift != round(shift)){
    stop("shift must be an integer.")
  }
  
  if(any(c(shape, rate, sigma) <= 0)){
    stop("shape, rate, and sigma must all be positive.")
  }
  
  if(!is.null(phiu)){
    if(phiu < 0 | phiu > 1){
      stop("phiu must be in [0,1].")
    }
  }
  
  if(is.null(phiu)){
    phiu <- 1-.pdiscgamma(u-1, shape=shape, rate=rate,
                         thresh=Inf, phiu = 0, shift=shift)
  }
  
  out <- rep(NA, length(x))
  
  if(sum(x>=u) != 0){
    out[x>=u] <- .ddiscgpd(x[x>=u], u, sigma, xi, phiu, log=log)
  }
  
  if(sum(x<u) != 0){
    out[x<u] <- .ddiscgamma(x[x<u], shape, rate, u, phiu, shift, log=log)
  }
  out
}

#' @importFrom VGAM dpareto
.eval_desponds <- function(t, Cminp, Cminq, alphap, alphaq){
  M <- 0.5*(dpareto(t, scale=Cminp, shape=alphap) +
              dpareto(t, scale=Cminq, shape=alphaq))
  
  one <- dpareto(t, scale=Cminp, shape=alphap)
  two <- dpareto(t, scale=Cminq, shape=alphaq)
  
  if(one == 0 & two == 0){
    out <- 0
  } else if(one == 0 & two != 0){
    out <- dpareto(t, scale=Cminq, shape=alphaq) *
      (log(dpareto(t, scale=Cminq, shape=alphaq))-log(M))
  } else if(one != 0 & two == 0){
    out <- dpareto(t, scale=Cminp, shape=alphap) *
      (log(dpareto(t, scale=Cminp, shape=alphap))-log(M))
  } else{
    out <- dpareto(t, scale=Cminp, shape=alphap) *
      (log(dpareto(t, scale=Cminp, shape=alphap))-log(M)) +
      dpareto(t, scale=Cminq, shape=alphaq) *
      (log(dpareto(t, scale=Cminq, shape=alphaq))-log(M))
  }
  out
}