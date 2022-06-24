#' Is object an ExpressionSet
#'
#' @param .exprs The object
#'
#' @return Logical if .exprs is an ExpressionSet
#' @export
#'
#' @examples
#' \dontrun{
#' is.ExpressionSet(affy::rma(Dilution))
#' }
is.ExpressionSet <- function(.exprs) {
  any(class(.exprs) %in% "ExpressionSet")
}
