# base R type check functions

.isNonEmptyDataFrame <- function(obj) {
    is.data.frame(obj) && sum(dim(obj)) > 0
}
assertthat::on.failure(.isNonEmptyDataFrame) <- function(call, env) {
    paste0(deparse(call$obj), " is not a non-empty `data.frame`")
}

.isListOfNonEmptyDataFrames <- function(obj) {
    is.list(obj) && all(sapply(obj, isNonEmptyDataFrame))
}
assertthat::on.failure(.isListOfNonEmptyDataFrames) <- function(call, env) {
    paste0(deparse(call$obj), " is not a list of non-empty `data.frame`s")
}

# bio objects

.is.seurat.object <- function(obj) inherits(obj, "Seurat")
assertthat::on.failure(.is.seurat.object) <- function(call, env) {
    paste0(deparse(call$obj), " is not a Seurat object")
}

.is.se.object <- function(obj) inherits(obj, "SummarizedExperiment")
assertthat::on.failure(.is.se.object) <- function(call, env) {
    paste0(deparse(call$obj), " is not a SummarizedExperiment object")
}

.is.seurat.or.se.object <- function(obj) {
    is.seurat.object(obj) || is.se.object(obj)
}
assertthat::on.failure(.is.seurat.or.se.object) <- function(call, env) {
    paste0(deparse(call$obj), " is not a Seurat or SummarizedExperiment object")
}

# general objects

.is.named.numeric <- function(obj) {
    is.numeric(obj) && !is.null(names(obj))
}
assertthat::on.failure(.is.named.numeric) <- function(call, env) {
    paste0(deparse(call$obj), " is not a named numeric vector")
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

.isIn <- function(x, table) {
    x %in% table
}

assertthat::on.failure(.isIn) <- function(call, env) {
    paste0(deparse(call$x), " is not in ", deparse(call$table))
}
