#' Load TCGA-OV Dataset
#'
#' Loads the TCGA-OV dataset from the "TCGA_OV.rda" file and returns it as a
#' SummarizedExperiment object.
#'
#' @return
#' A SummarizedExperiment object containing raw RNA-Seq counts from the TCGA-OV dataset
#' with associated metadata columns.
#'
#' @usage
#' load_TCGA_OV()
#'
#' @examples
#' # Load the TCGA-OV dataset
#' my_dataset <- load_TCGA_OV()
#'
#' # Access metadata columns
#' metadata(my_dataset)$AGE
#' metadata(my_dataset)$TUMOR_GRADE
#'
#' @export
load_TCGA_OV <- function() {
  load("TCGA_OV.rda")
  if (exists("TCGA_OV", inherits = FALSE)) {
    return(TCGA_OV)
  } else {
    stop("TCGA_OV not found in the loaded data.")
  }
}
