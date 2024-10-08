% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate.R
\name{annotation_object}
\alias{annotation_object}
\title{Return Probeset Annotation Object for ExpressionSet}
\usage{
annotation_object(.exprs)
}
\arguments{
\item{.exprs}{The expression set that has a corresponding (affy) annotation.}
}
\value{
The annotation.db associated with the ExpressionSet annotation.
}
\description{
Get the right annotation (possibly installing it) for an expressionset.
}
\details{
This is a fragile function. It takes either a \link[Biobase:class.ExpressionSet]{Biobase::ExpressionSet}
or \link[AnnotationDbi:AnnotationDb-class]{AnnotationDbi::ChipDb} and returns the \link[AnnotationDbi:AnnotationDb-class]{AnnotationDbi::ChipDb}.
That way, a user can provide the expression set and still get the
desired annotation. This could be object
oriented, but it didn't seem worth the trouble/cost of polluting the namespace.

If the parameter is an ExpressionSet, the function looks at the
annotation slot of the \link[Biobase:class.ExpressionSet]{Biobase::ExpressionSet} and tries to install the
corresponding annotation package for it. It then returns the corresponding
annotation object (db) file to be used in annotation lookups.

If the parameter is a \link[AnnotationDbi:AnnotationDb-class]{AnnotationDbi::ChipDb}, it just returns it.
}
\examples{
\dontrun{
library(affydata)
data(Dilution)
# Get the first 5 probesets of the annotation.
annotation_object(affy::rma(Dilution)) |>
   AnnotationDbi::keys() |>
   magrittr::extract(1:5)
}
}
