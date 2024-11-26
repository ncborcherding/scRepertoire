# base R type check functions

isNonEmptyDataFrame <- function(obj) {
    is.data.frame(obj) && sum(dim(obj)) > 0
}
assertthat::on_failure(isNonEmptyDataFrame) <- function(call, env) {
    paste0(deparse(call$obj), " is not a non-empty `data.frame`")
}

isListOfNonEmptyDataFrames <- function(obj) {
    is.list(obj) && all(sapply(obj, isNonEmptyDataFrame))
}
assertthat::on_failure(isListOfNonEmptyDataFrames) <- function(call, env) {
    paste0(deparse(call$obj), " is not a list of non-empty `data.frame`s")
}

# bio objects

is_seurat_object <- function(obj) inherits(obj, "Seurat")
assertthat::on_failure(is_seurat_object) <- function(call, env) {
    paste0(deparse(call$obj), " is not a Seurat object")
}

is_se_object <- function(obj) inherits(obj, "SummarizedExperiment")
assertthat::on_failure(is_se_object) <- function(call, env) {
    paste0(deparse(call$obj), " is not a SummarizedExperiment object")
}

is_seurat_or_se_object <- function(obj) {
    is_seurat_object(obj) || is_se_object(obj)
}
assertthat::on_failure(is_seurat_or_se_object) <- function(call, env) {
    paste0(deparse(call$obj), " is not a Seurat or SummarizedExperiment object")
}

# general objects

is_named_numeric <- function(obj) {
    is.numeric(obj) && !is.null(names(obj))
}
assertthat::on_failure(is_named_numeric) <- function(call, env) {
    paste0(deparse(call$obj), " is not a named numeric vector")
}

# functions

isIn <- function(x, table) {
    x %in% table
}

assertthat::on_failure(isIn) <- function(call, env) {
    paste0(deparse(call$x), " is not in ", deparse(call$table))
}
