#' Subset an ExpressionSet by a list of identifiers
#'
#' @description Subset an ExpressionSet by providing a list of probe/gene identifiers.
#'
#' @details This function subsets an ExpressionSet by probe (row). This can be done
#' directly from the [Biobase::ExpressionSet] object, however there are a few additional
#' features of this function.
#'
#' First, you can use other fields (like SYMBOL) which are present in annotation to do
#' a subsetting. It handles the mappings back and forth to probesets.
#'
#' Second, it tells you what it's doing. Not all of your identifiers may be present (this
#' happens in gene signatures, etc). So it will tell you that it's happening (with `verbose=TRUE`).
#'
#' @param .exprs An ExpressionSet to subset
#' @param subset A list of probeset identifiers, or other identifier (e.g., gene symbol)
#' @param by A selector of the mapping (default is PROBESET, choices include SYMBOL).
#' @param verbose Should statistics be reported on what is subset.
#' @param ... Any additional arguments
#'
#' @return A new ExpressionSet with fewer rows, representing a subset.
#' @export
#'
#' @examples
#' \dontrun{
#' # Subset by probesets
#' subset_by(exprs, subset=c("1007_s_at","1008_s_at"))
#' # Subset by gene
#' subset_by(exprs, subset=c("STAT1","EGFR","TP53"), by = "SYMBOL")
#' }
subset_by <- function(.exprs, subset = NULL, by = c("PROBESET","SYMBOL"), verbose = TRUE) {
  by = match.arg(by, c("PROBESET","SYMBOL"))

  stopifnot(is.ExpressionSet(.exprs))

  names(subset)<-subset
  ps <- Biobase::featureNames(.exprs)

  # Based on how to subset, the selection of keepers, missing and is different.
  if ( by == "PROBESET") {
    ps_keepers <- intersect(ps, subset)
    subset_keepers <- ps_keepers
    subset_missing <- setdiff(subset, ps)
  } else if ( by == "SYMBOL" ) {
    s <- annotate_with_symbol(.exprs, use = "everything")
    matches <- purrr::map_lgl(s, ~any(.x %in% subset))
    ps_keepers <- names(which(matches))
    m <- subset %in% unlist(s)
    subset_keepers <- subset[m]
    subset_missing <- subset[!m]
  }

  if (!by %in% c("PROBESET","SYMBOL"))
    stop("NO OTHER BY OPTION IMPLEMENTED YET!")

  cli::cli_inform(c(
    "subset_by: subsetting ExpressionSet by {by} ({length(subset)} identifiers).",
    "i" = glue::glue("{length(subset_missing)} identifiers not in ExpressionSet."),
    "i" = glue::glue("{length(ps)-length(ps_keepers)} probesets removed from ExpressionSet."),
    "i" = glue::glue("{length(subset_keepers)} identifiers remaining in ExpressionSet ({length(ps_keepers)} probesets).")
  ))

  .exprs[ps_keepers,]

}

#' @describeIn subset_by Subset by gene symbol
#' @export
subset_by_symbol <- function(...) {
  subset_by(by="SYMBOL",...)
}
