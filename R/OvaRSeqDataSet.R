#' OvaRSeqDataSet class
#'
#' \code{OvaRSeqDataSet} is a class inherited from
#' \code{\link{SummarizedExperiment}}.
#' It is used to store the count matrix, colData, and design formula
#' in differential expression analysis.
#' @references Martin Morgan, Valerie Obenchain, Jim Hester and
#' Hervé Pagès (2018). SummarizedExperiment: SummarizedExperiment container.
#' R package version 1.12.0.
#' @export OvaRSeqDataSet
setClass("OvaRSeqDataSet",
         contains="SummarizedExperiment",
         representation=representation(
           design = "ANY"
         )
)

#' OvaRSeq constructor
#' @param countData a matrix or data frame contains gene count
#' @param colData a \code{DataFrame} or \code{data.frame}
#' @param ... optional arguments passed to \code{SummarizedExperiment}
#' @return a OvaRSeq object
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
OvaRSeqDataSet <- function(countData, colData){

  # coerce count data to matrix
  count_data <- as.matrix(count_data)

  # coerce col data to DataFrame
  col_data <- DataFrame(col_data)

  # check that count data and col data have same number of columns
  stopifnot(ncol(count_data) == nrow(col_data))

  # create SummarizedExperiment object
  se <- SummarizedExperiment(assays = list(counts = count_data), colData = col_data, ...)
  object = new("OvaRSeqDataSet", se)
  object
}
