
#' Translate ExpressionSet to alternate dimensions.
#'
#' @description An ExpressionSet can be morphed into a new dimensionality, but only carefully.
#'
#' @details An ExpressionSet needs to maintain the correct dimensions for the linked data
#' within it. This function attempts to clone an ExpressionSet then alter one of the dimensions
#' based on input data. All other data is taken from the original copy.
#'
#' @param src Original ExpressionSet to clone
#' @param new_exprs A new matrix of expression data (with either fewer rows or columns)
#' @param pdata Potentially new phenoData to use
#' @param fdata Potentially new featureData to use
#' @param annotation Potentiall new annotation
#'
#' @return A new expression set
#' @export
#'
#' @examples
#' \dontrun{
#' translate_ExpressionSet(exprset, matrix(c(1,2,3)))
#' }
translate_ExpressionSet <- function(src, new_exprs, pdata = NULL, fdata = NULL, annotation = NULL) {

    x <- Biobase::ExpressionSet(
    assayData = new_exprs,
    phenoData = (
      if (!is.null(pdata))
        Biobase::AnnotatedDataFrame(pdata)
      else if (ncol(src) == ncol(new_exprs))
        Biobase::phenoData(src)
      else
        Biobase::annotatedDataFrameFrom(new_exprs, byrow = FALSE)
    ),
    featureData = (
      if (!is.null(fdata))
        Biobase::AnnotatedDataFrame(fdata)
      else if (nrow(src) == nrow(new_exprs))
        Biobase::featureData(src)
      else
        Biobase::annotatedDataFrameFrom(new_exprs, byrow = TRUE)
    ),
    experimentData = Biobase::experimentData(src),
    protocolData = Biobase::protocolData(src),
    annotation = annotation
  )
}
