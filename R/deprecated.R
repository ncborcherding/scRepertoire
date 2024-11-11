#' ***DEPRECATED*** Take the meta data in seurat/SCE and place it into a list
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' Allows users to perform more fundamental measures of clonotype analysis
#' using the meta data from the seurat or SCE object. For Seurat objects the
#' active identity is automatically added as "cluster". Remaining grouping
#' parameters or SCE or Seurat objects must appear in the meta data.
#'
#' This function is deprecated as of version 2 due to the confusion it caused
#' to many users. Users are encouraged to remain with the abstraction barrier
#' of combined single cell objects and the outputs of [combineTCR()] and
#' [combineBCR()] for all functions.
#'
#' We discourage the use of this function, but if you have to use it, set the
#' `force` argument to `TRUE`.
#'
#' @param sc output of [combineExpression()].
#' @param ... previously the `group` or `split.by` argument, indicating the
#' column header to group the new list by. This should strictly be one argument
#' and is an ellipsis for backwards compatibility. Everything after the first
#' argument is ignored.
#' @param force logical. If not `TRUE` (default), a deprecation error will
#' be thrown. Otherwise the function will run but not guaranteed to be stable.
#'
#' @export
#' @return list derived from the meta data of single-cell object with
#' elements divided by the group parameter
#' @keywords internal
#'
expression2List <- function(sc, ..., force = FALSE) {

    if (!force) {
        lifecycle::deprecate_stop(
            "2.0.0", "expression2List()",
            details = "If you MUST use it, see `?expression2List`"
        )
    }

    lifecycle::deprecate_soft("2.0.0", "expression2List()")
    .expression2List(sc, list(...)[[1]])
}
