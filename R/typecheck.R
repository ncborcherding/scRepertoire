# scRepertoire objects

isCombineContigsOutput <- function(obj) {
    is.list(obj) && all(sapply(obj, is.data.frame))
}
assertthat::on_failure(isCombineContigsOutput) <- function(call, env) {
    paste0(deparse(call$obj), " is not an output of combineTCR or combineBCR")
}

isListOfTwoCombineContigsOutputs <- function(obj) {
    is.list(obj) && length(obj) == 2 && all(sapply(obj, isCombineContigsOutput))
}
assertthat::on_failure(isListOfTwoCombineContigsOutputs) <- function(call, env) {
    paste0(
        deparse(call$obj),
        " is not a list of two outputs of combineTCR and combineBCR"
    )
}

isAnyValidProductOfCombineContigs <- function(obj) {
    isCombineContigsOutput(obj) || isListOfTwoCombineContigsOutputs(obj)
}
assertthat::on_failure(isAnyValidProductOfCombineContigs) <- function(call, env) {
    paste0(
        deparse(call$obj),
        " is not a valid output of combineTCR or combineBCR, nor a list of them"
    )
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
