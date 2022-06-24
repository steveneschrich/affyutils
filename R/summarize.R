#' Summarize ExpressionSet rows to reduced form
#'
#' @description Summarize an ExpressionSet to a reduced form by grouping
#' on a feature annotation variable.
#'
#' @details We often need to summarize an ExpressionSet by Gene Symbol or other
#' annotation field. There are two parts to having to do this work. First, you
#' have to map/group the probesets by this annotation field. Then you have to
#' pick a summarization method. This function does those two things and
#' returns an ExpressionSet with the reduced information, copying where
#' possible the attributes of the original ExpressionSet.
#'
#' @param .exprs An ExpressionSet to summarize
#' @param by The annotation field to summarize by (e.g., SYMBOL)
#' @param method The method to use for summarization. Can include median, max, mean, min, tukey.biweight.
#' @param ... Any additional parameters
#'
#' @return A summarized ExpressionSet containing phenodata, etc from the original ExpressionSet
#' but with summarize rows.
#'
#' @export
#' @importFrom rlang .data
#' @examples
#' \dontrun{
#' summarize_by(exprset, by = "SYMBOL", method = "max")
#' }
summarize_by <- function(.exprs, by = "SYMBOL", method = c("median","max","mean","min", "tukey.biweight")) {
  stopifnot(is.ExpressionSet(.exprs))

  method <- match.arg(method, c("median","max","mean","min", "tukey.biweight"))
  s <- purrr::discard(annotate_with(.exprs, with = by), is.na)

  # Using some dplyr to get indices of same symbol
  mapping <- tibble::enframe(s, name=NULL) |>
    dplyr::mutate(n=dplyr::row_number()) |>
    dplyr::group_by(.data$value) |>
    dplyr::summarize(indices = list(.data$n)) |>
    tibble::deframe()

  summarizer <- switch(method,
                       median = .summarize_median,
                       mean = .summarize_mean,
                       max = .summarize_max,
                       min = .summarize_min,
                       tukey.biweight = .summarize_tukey_biweight
  )
  new_exprs <- purrr::map_dfr(mapping, summarizer, Biobase::exprs(.exprs), .id = "rn") |>
    tibble::column_to_rownames("rn") |>
    as.matrix()

  # TODO: Annotate this expression set with gene symbol information somehow.
  .y <- translate_ExpressionSet(.exprs, new_exprs = new_exprs, annotation=by)


  .y
}

#' @describeIn summarize_by Summarize by gene symbol.
#' @export
summarize_by_symbol <- function(...) {
  summarize_by(by="SYMBOL",...)
}

#' @describeIn summarize_by Summarize by entrez gene ID.
#' @export
summarize_by_entrezid <- function(...) {
  summarize_by(by = "ENTREZID", ...)
}

#' @describeIn summarize_by Summarize by Ensemble ID.
#' @export
summarize_by_ensembl <- function(...) {
  summarize_by(by = "ENSEMBL", ...)
}

#' Summarize a group of probesets
#'
#' @name summarize
#'
#' @description Given a set of samples and a set of probesets, summarize to one measurement
#'  per sample.
#'
#' @details Summarizing data from multiple probesets to a single value can be tricky since not
#' all probesets are equally good.
#'
#' @param i The rows of `x` to use in summarizing
#' @param x The matrix of expression (rows=genes, columns=samples)
#'
#' @return A summarized matrix with same number of columns, but one row.
#'
NULL

#' @describeIn summarize Use the max value
#' @export
.summarize_max <- function(i, x) {
  apply(x[i,,drop=FALSE], 2, max)
}

#' @describeIn summarize Use the min value
#' @export
.summarize_min <- function(i, x) {
  apply(x[i,,drop=FALSE], 2, min)
}

#' @describeIn summarize Use the median value
#' @export
.summarize_median <- function(i, x) {
  apply(x[i,,drop=FALSE], 2, stats::median)

}

#' @describeIn summarize Use the mean value
#' @export
.summarize_mean <- function(i, x) {
  colMeans(x[i,,drop=F])
}

.summarize_tukey_biweight <- function(i, x) {
  affy::tukey.biweight(x[i,,drop=F])
}


