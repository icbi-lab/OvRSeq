#' Build SummarizedExperiment Object for OvrSeq
#'
#' This function creates a SummarizedExperiment object for OvrSeq analysis by combining count data and column metadata.
#'
#' @param count_data_df A data frame containing count data (log2-transformed TPM or similar).
#' The rows should represent genes, and columns should represent samples, with row names set as gene names.
#' @param col_data_df A data frame containing column metadata.
#' Each row should correspond to a sample, and columns should represent sample attributes. Row names should be set as sample IDs.
#'
#' @return A SummarizedExperiment object with the count data and column metadata integrated.
#'
#' @details
#' The function ensures that only samples present in both count data and column metadata are used to construct the SummarizedExperiment object. It checks for the presence of necessary row names in both data frames and then performs the intersection. The count data is transposed to match the format expected by the SummarizedExperiment constructor.
#'
#' @examples
#' # Assuming count_data and col_data are pre-loaded data frames
#' se <- build_OvrSeq_SE(count_data, col_data)
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
OvaRSeqDataSet <- function(count_data_df, col_data_df) {
  # Check if row names are present in count_data_df
  if (is.null(rownames(count_data_df))) {
    stop("Row names (gene names) are required in count data frame.")
  }

  # Check if row names are present in col_data_df
  if (is.null(rownames(col_data_df))) {
    stop("Row names (sample IDs) are required in column data frame.")
  }

  # Intersect gene names and sample IDs between countData and colData
  common_samples <- intersect(rownames(col_data_df), colnames(count_data_df))

  # Subset the data frames to include only common genes and samples
  count_data_sub <- count_data_df[, common_samples]
  col_data_sub <- col_data_df[common_samples, ]

  # Create SummarizedExperiment object
  se <- SummarizedExperiment(assays = list(counts = as.matrix(count_data_sub)),
                             colData = as.data.frame(col_data_sub))

  # Return the SummarizedExperiment object
  return(se)
}
