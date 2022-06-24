#' Convert a one column data frame with rownames to a list
#'
#' @details The [tibble::deframe()] function is awesome, but since rownames are always
#' special, it doesn't quite work as I'd like. This function will transform a single
#' column data frame (with rownames) into a named list. It appears that [base::drop()]
#' does not honor row names so this seems the best option.
#'
#' @note The data frame column name is lost.
#'
#' @param .x A data frame
#'
#' @return A list with rownames as the list names, and the values as the data frame column.
#' @export
#'
#' @importFrom rlang .data
#' @examples
#' \dontrun{
#' # Iris has numbers for rownames, so this will have number=value
#' deframe_with_rownames(data.frame(iris[,1]))
#' }
deframe_with_rownames<-function(.x) {
  stopifnot(ncol(.x)==1)
  .x |>
    tibble::rownames_to_column(var = "rn") |>
    dplyr::select(.data$rn, dplyr::everything()) |>
    tibble::deframe()
}


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
