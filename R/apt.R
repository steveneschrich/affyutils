#' Load APT data into ExpressionSet
#'
#' @description Parse and load APT output into a [[Biobase::ExpressionSet()]] object.
#'
#' @details APT produces two output files of interest: the report and summary files. The
#' report file consists of a series of lines of parameter settings (prefixed with `#%`)
#' and rows of statistics about the CEL file that was processed. The summary file starts
#' with the parameter settings (same as report file) and includes the gene expression estimates.
#'
#' This function will load both files (using `basename` as the file stem for both the
#' summary and report files). The information will be stored in a [[Biobase::ExpressionSet()]]
#' object. There are a couple of points about the object:
#'
#' - `assayData` will contain the gene expression values
#' - `annotation` will contain the chip name from the metadata
#' - `pData` will contain the statistics about the chip (average intensity, etc)
#' - `preproc(experimentData())` contains the algorithm parameters used in the APT process
#'
#' @param basename The file stem to use for loading the summary and report txt files
#'
#' @return A [[Biobase::ExpressionSet()]] object representing the contents of the file(s).
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' # Load rma.summary.txt and rma.report.txt
#' load_apt_rma(basename = "rma")
#' }
load_apt_rma <- function(basename="rma") {
  summary_file <- paste0(basename, ".summary.txt")

  report <- load_apt_rma_report(paste0(basename, ".report.txt"))

  summary <- load_apt_rma_summary(paste0(basename, ".summary.txt"))

  stopifnot(report$statistics$cel_files == colnames(summary))

  Biobase::pData(summary)<- report$statistics

  summary
}



#' Load the APT RMA summary file
#'
#' @param f The file containing RMA summary information
#'
#' @return A [[Biobase::ExpressionSet()]] with the expression data (assayData) and experimentInfo.
#'
#' @examples
load_apt_rma_summary <- function(f) {

  # The header is proportional to the number of samples that were processed. There is no
  # good way of knowing how big the header is, since we don't know how many samples. Even
  # the file size is not accurate, since the chip may have different # of probes. So we
  # just try multiple lengths until finding enough. This file can be really big, so this
  # is presumably still a better option.
  all_header <- FALSE
  n <- 200
  while (!all_header) {
    hdr <- readLines(f, n = n<-n*2)
    all_header <- index_of_last_apt_comment(hdr) != length(hdr)
  }



  exprs <- readr::read_tsv(f,
                           skip = index_of_last_apt_comment(hdr),
                           col_types = readr::cols("probeset_id"="c",.default = "d")
  ) |>
    tibble::column_to_rownames("probeset_id") |>
    as.matrix()

  params <- parse_apt_header(hdr)

  Biobase::ExpressionSet(
    assayData = exprs,
    experimentData = Biobase::MIAME(
      preprocessing = as.list(params)
    ),
    annotation = extract_apt_platform_string(hdr)
  )
}

#' Extract header lines from a vector of strings
#'
#' @description Extract header lines from a vector of lines from a file.
#'
#' @details APT headers consist of the prefix `#%`. This routine finds these
#' lines, strips the prefix and returns the result as a list.
#'
#' @param hdr A vector of strings to examine
#'
#' @return A list of header name=value pairs
#'
#' @importFrom rlang .data
parse_apt_header <- function(hdr) {

  algorithm_parameters <- hdr |>
    tibble::enframe(name = NULL) |>
    dplyr::filter(stringr::str_starts(.data$value,"#%")) |>
    dplyr::mutate(value = stringr::str_remove(.data$value, "^#%")) |>
    dplyr::pull("value")
}


#' Find row number of last comment in vector of strings
#'
#' @param s A vector of strings
#'
#' @return The row number of the last comment
#'
#' @examples
index_of_last_apt_comment <- function(s) {
  max(which(stringr::str_starts(s,"#%")))
}


#' Parse and return the report file from APT
#'
#' @param f A file with the report content.
#'
#' @return A list with algorithm parameters (`parameters`) and normalization statistics (`statistics`).
#'
#' @importFrom rlang .data
load_apt_rma_report <- function(f) {
  content <- readLines(f)

  algorithm_parameters <- parse_apt_header(content)
  stats <- content |>
    tibble::enframe(name = NULL) |>
    dplyr::filter(!stringr::str_starts(.data$value, "#%")) |>
    dplyr::pull("value") |>
    I() |>
    readr::read_tsv( col_types = readr::cols("cel_files"="c",.default = "d"))

  list(parameters = algorithm_parameters, statistics = stats)
}

#' Extract platform info from APT header
#'
#' @description Given a vector of strings, extract the platform or chip type.
#'
#' @details In the APT output, the parameter
#' ```
#' affymetrix-algorithm-param-apt-opt-chip-type
#' ```
#' describes the Affymetrix platform (chip) used in the analysis. This routine will
#' find this entry in the vector and extract it.
#'
#' @param s A vector of strings representing APT parameters
#'
#' @return A string representing the extracted string type
#'
#'
#' @examples
#' \dontrun{
#' extract_apt_platform_string(c("#%affymetrix-algorithm-param-apt-opt-chip-type=Clariom-D"))
#' }
extract_apt_platform_string <- function(s) {

  platform_string <- stringr::str_subset(s, "affymetrix-algorithm-param-apt-opt-chip-type=") |>
    unique() |>
    stringr::str_match("chip-type=(.*)") |>
    magrittr::extract(2)

  platform_string
}

