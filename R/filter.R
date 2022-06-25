#' Filter genes with range less than threshold
#'
#' @param .exprs ExpressionSet to use for filtering.
#' @param threshold Range to use as minimum (default 1)
#'
#' @return A new ExpressionSet with low range genes removed.
#' @export
#'
#' @examples
filter_low_range<-function(.exprs, threshold=1) {
  stopifnot(class(.exprs) %in% "ExpressionSet")

  ranges<-Biobase::rowMax(.exprs) - Biobase::rowMin(.exprs)

  .exprs[which(ranges>threshold),]
}


#' Filter out genes (features) with no samples above noise (low expression)
#'
#' @param .exprs the ExpressionSet to use
#' @param threshold  a minimum acceptable threshold of gene expression (default 5)
#'
#' @return A new ExpressionSet reduced to genes passing the criteria.
#' @export
#'
#' @examples
filter_low_expression<-function(.exprs, threshold=5) {
  .exprs[which(Biobase::rowMax(.exprs)>threshold),]
}

#' Filter genes on number of samples above a target threshold
#'
#' @param .exprs ExpressionSet
#' @param threshold A threshold to consider low
#' @param n The number of samples exceeding threshold
#'
#' @note If n is fractional, it is treated as a percentage.
#'
#' @return ExpressionSet reduced to rows with more than `n` columns larger than `threshold`.
#' @export
#'
#' @examples
#' \dontrun{
#' # Filter probesets with more than 20% above expression level of 5.
#' filter_n_above_threshold(exprset, threshold = 5, n = 0.2)
#'
#' # Filter probesets with more than 10 samples above expression level of 10.
#' filter_n_above_threshold(exprset, threshold = 10, n = 10)
#' }
filter_n_above_threshold <- function(.exprs, threshold = 5, n=0) {
  if ( n > 0 && n < 1 )
    n <- ceiling(n*ncol(.exprs)[[1]])

  num_above <- apply(Biobase::exprs(.exprs), 1, function(x){length(which(x>threshold))})
  .exprs[which(num_above >= n),]
}

#' Filter out genes (features) with low variability
#'
#' @param .exprs the ExpressionSet to use
#' @param threshold  a minimum acceptable threshold of gene expression variance
#'
#' @return A new ExpressionSet reduced to genes passing the criteria.
#' @export
#'
#' @examples
#' \dontrun{
#' filter_by_variance(exprset, threshold = 0.5)
#' }
filter_by_variance <- function(.exprs, threshold = 0 ) {
  vars <- apply(Biobase::exprs(.exprs), 1, stats::var)
  .exprs[which(vars >= threshold),]
}
