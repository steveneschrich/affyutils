# NB: I don't like prannot, that needs to change.



#' Annotate Probesets of ExpressionSet
#'
#' @description Add [Biobase::featureData] to a [Biobase::ExpressionSet] consisting of probeset annotation.
#'
#' @details This function is inspired by affycoretools. If the standard approach when building an
#' ExpressionSet (from affy data) is to include annotation then downstream work is much easier. This
#' is that approach.
#'
#' @note There is a great discussion in affycoretools about the problem of one-to-many and many-to-many
#' mapping of annotation. This package uses the first entry approach, which is not ideal but is practical.
#'
#' @param .exprs A [Biobase::ExpressionSet] to annotate
#'
#' @return A [Biobase::ExpressionSet] with [Biobase::featureData] added incorporating existing chip annotation.
#' @export
#'
#' @examples
#' \dontrun{
#' annotate_exprset(exprs)
#' }
annotate_exprset <- function(.exprs) {

  stopifnot(is.ExpressionSet(.exprs))
  # This would be great if it is a named list, with many or one
  # so you could put in c("ACCNUM"="list", "SYMBOL"="first") and it
  # would do the appropriate thing per field. Interesting.
  #a <- get_annotation(.exprs)

  #AnnotationDbi::mapIds(a, keys = Biobase::featureNames(.exprs), column = "", keytype = "PROBEID")
}

#' Return Available Annotation Fields
#'
#' @param .exprs A [Biobase::ExpressionSet] or [AnnotationDbi::ChipDb-class]
#'
#' @return The list of columns in the annotation data.
#' @export
#'
#' @examples
annotation_fields <- function(.exprs) {
  AnnotationDbi::columns(annotation_object(.exprs))
}

#' Generate Map of Probeset Annotation
#'
#' @details This is a composite function. In its simplest form, it
#' will provide a map of probesets to an annotation fields (see
#' [annotation_fields()] for fields). It does a simple approach defaulted
#' to in [AnnotationDbi::mapIds] and selects the first annotation. This map
#' will be a named vector with the name being the probeset and the value
#' being the mapped field.
#'
#' If multiple fields are provided, then this process is repeated for all
#' fields and the result is a data frame with rownames as probeset and each
#' field being a separate column. Note that the "take-the-first" approach
#' still holds for each column, which could provide some surprising results.
#'
#' @param .exprs An [Biobase::ExpressionSet] or [AnnotationDbi::ChipDb-class] object
#' @param with The fields to map from probesets
#' @param use The selection process when more than one element matches. See [selectors()]
#' for details.
#' @param ... Any further arguments to the function
#'
#' @return A named vector or data frame with names/rownames as probesets
#'    and the mapped data as values.
#'
#' @export
#'
#' @examples
annotate_with <- function(.exprs, with=NULL, use=c("first","last","everything")) {
  a <- annotation_object(.exprs)
  use <- match.arg(use,c("first","last","everything"))
  if (use == "everything")
    f <- .map_fields_12many(a, fields = with)
  else
    f <- .map_fields_121(a, fields = with)

  # IMPORTANT: we have to put annotation in the exprs probeset order!
  if ( is.null(dim(f)) )
    f[Biobase::featureNames(.exprs)]
  else
    f[Biobase::featureNames(.exprs),]
}

#' @describeIn annotate_with Annotate with gene symbol
#' @export
annotate_with_symbol <- function(...) {
  annotate_with(..., with = "SYMBOL")
}

#' @describeIn annotate_with Annotate with Ensembl ID
#' @export
annotate_with_ensembl <- function(...) {
  annotate_with(..., with = "ENSEMBL")
}

#' @describeIn annotate_with Annotate with entrez gene ID
#' @export
annotate_with_entrezid <- function(...) {
  annotate_with(..., with = "ENTREZID")
}


#' Select an annotation element
#'
#' @description Given a list of annotation elements, select the most
#'   representative one.
#'
#' @details In one-to-many mappings, there are many ways to pick a representative
#' element. It could be the first, last, etc. These routines are selectors that
#' implement these various logical choices.
#' @param .x A list of character vectors to select from.
#' @name selectors
#'
NULL

#' @describeIn selectors Select the first element.
.select_first_element <- function(.x) {.x[[1]]}
#' @describeIn selectors Select the last element.
.select_last_element <- function(.x) {.x[[length(.x)]]}
#' @describeIn selectors Select the element with the longest string length.
.select_longest_element <- function(.x) {
  z <- purrr::map_int(.x, stringr::str_length)
  .x[[which(z==max(z))[1]]]
}


#' Generate Map of Annotation
#'
#' @description Maps annotation data from probeset to other columns but
#' chooses only one match.
#'
#' @details This is an annotation function that maps from probesets to
#' other fields. It selects only one mapping, even if multiple mappings
#' exist. The selector parameter determines how that map is chosen.
#'
#' @param a An [AnnotationDbi::ChipDb-class] object.
#' @param fields The fields to extract
#' @param selector A selector when multiple matches are found
#'
#' @return A list or data frame of mapped data with names being the probesets.
#' @export
#'
#' @examples
.map_fields_121 <- function(a,fields = NULL, selector = .select_first_element ) {

  names(fields)<-fields
  res <- purrr::map_dfr(fields, ~AnnotationDbi::mapIds(a, keys = AnnotationDbi::keys(a),
                                                column = ., keytype = "PROBEID",
                                                multiVals = selector),
                 .id = "FIELD") |>
    tibble::column_to_rownames("FIELD") |>
    t() |>
    data.frame()

  # Special case one column data, convert to list.
  if ( ncol(res) == 1)
    res <- saeutils::deframe_with_rownames(res)

  res
}




#' Generate One-to-Many Map of Probeset Annotation
#'
#' @description Provides a mapping of annotation as embedded lists.
#'
#' @details This is a somewhat useful function, but be careful what you wish for. This
#' function provides the one-to-many mapping of annotation via a data frame
#' full of lists. That is, for a given probe and a given field it will be
#' a list of values. For something like accession numbers this could be a list
#' of 50 accessions. Embedded within a single position in the data frame of
#' probeset x field. Still, if you want that level of detail then this is the
#' object for you.
#'
#'
#' @param .exprs An ExpressionSet or Annotation Object
#' @param fields The fields to extract
#' @return A data frame of probes by fields with many values embedded as lists.
#' @export
#'
#' @examples
.map_fields_12many <- function(.exprs, fields = NULL) {

  a <- annotation_object(.exprs)

  names(fields)<-fields
  res <- purrr::map(fields, ~AnnotationDbi::mapIds(a, keys = AnnotationDbi::keys(a),
                                                column = ., keytype = "PROBEID",
                                                multiVals = "list")
  ) |>
    purrr::map_dfr(~purrr::map(., ~list(.)), .id = "FIELD") |>
    tibble::column_to_rownames("FIELD") |>
    t() |>
    data.frame()

  # Special case one column data, convert to list.
  if ( ncol(res) == 1)
    res <- saeutils::deframe_with_rownames(res)

  res
}




#' Return Probeset Annotation Object for ExpressionSet
#'
#' @description Get the right annotation (possibly installing it) for an expressionset.
#'
#' @details This is a fragile function. It takes either a [Biobase::ExpressionSet]
#' or [AnnotationDbi::ChipDb-class] and returns the [AnnotationDbi::ChipDb-class].
#' That way, a user can provide the expression set and still get the
#' desired annotation. This could be object
#' oriented, but it didn't seem worth the trouble/cost of polluting the namespace.
#'
#' If the parameter is an ExpressionSet, the function looks at the
#' annotation slot of the [Biobase::ExpressionSet] and tries to install the
#' corresponding annotation package for it. It then returns the corresponding
#' annotation object (db) file to be used in annotation lookups.
#'
#' If the parameter is a [AnnotationDbi::ChipDb-class], it just returns it.
#'
#' @param .exprs The expression set that has a corresponding (affy) annotation.
#'
#' @return The annotation.db associated with the ExpressionSet annotation.
#' @export
#'
#' @examples
#' \dontrun{
#' library(affydata)
#' data(Dilution)
#' # Get the first 5 probesets of the annotation.
#' annotation_object(affy::rma(Dilution)) |>
#'    AnnotationDbi::keys() |>
#'    magrittr::extract(1:5)
#' }
annotation_object<-function(.exprs) {

  if (is.ExpressionSet(.exprs)) {
    annotationName<-Biobase::annotation(.exprs)
    annotationPackage<-paste0(annotationName,".db")
    if (!base::requireNamespace("BiocManager", quietly = TRUE))
      utils::install.packages("BiocManager")
    if (!base::requireNamespace(annotationPackage, quietly=TRUE))
      BiocManager::install(annotationPackage)

    # Return the annotation object (for select, etc). It is the name of the package.
    a <- base::getExportedValue(annotationPackage, annotationPackage)
  } else if (class(.exprs) %in% "ChipDb")
    a <- .exprs
  else
    stop("Cannot handle object of class ", class(.exprs))

  a
}



