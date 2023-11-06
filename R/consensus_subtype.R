
#' Obtain Consensus Subtypes of High-Grade Serous Ovarian Cancer
#'
#' This function gets the consensus molecular subtypes of high-grade serous (HGS) ovarian cancer using the `consensusOV` package.
#' The function implements a consensus classifier which consolidates and improves on the robustness of proposed subtype classifiers.
#'
#' @param se A Summarized Experiments of gene expression values, with gene symbol or entrez id as identifiers.
#' @param ids_type A character string indicating the type of IDs used in the 'data' row names. Can be either "symbol" or "entrez".
#' @param concordant.tumors.only Logical. If TRUE, only tumors that are concordantly classified across all datasets are included in the final classification. Defaults to TRUE.
#' @param remove.using.cutoff Logical. If TRUE, tumors with poor classification confidence are removed using the optimal cutoff. Defaults to FALSE.
#' @param percentage.dataset.removed Numeric value indicating the percentage of the dataset to be removed based on classification confidence. Defaults to 0.75.
#' @return Summarized Experiments with Tumor_Molecular_Subtypes of consensus molecular subtypes for the samples, in colData.
#' @importFrom consensusOV get.consensus.subtypes
#' @examples
#' load_TCGA_OV
#' TCGA_OV <- get_consensus_ov_subtypes(TCGA_OV, ids_type = "symbol")
#' @export
#' @references
#' Chen G, Kannan L, Geistlinger L, Kofia V, Safikhani Z, Gendoo D, Parmigiani G, Birrer M, Haibe-Kains B, Waldron L (2018).
#' "Consensus on molecular subtypes of high-grade serous ovarian carcinoma."
#' Clinical Cancer Research, 24, 4990. doi:10.1158/1078-0432.CCR-18-0784.
get_consensus_ov_subtypes <- function(se, ids_type = "symbol",
                                      concordant.tumors.only = TRUE, remove.using.cutoff = FALSE,
                                      percentage.dataset.removed = 0.75) {

  if (ids_type == "symbol") {
    # Convert gene symbols to Entrez IDs if needed
    ids <- update_se_with_entrez_ids(se)
  } else if (ids_type == "entrez") {
    ids <- rownames(data)
  } else {
    stop("Invalid 'ids_type'. Choose either 'symbol' or 'entrez'.")
  }


  # Obtain consensus subtypes using the get.consensus.subtypes function from consensusOV
  subtypes <- get.consensus.subtypes(assay(se), ids,
                                     concordant.tumors.only = concordant.tumors.only,
                                     remove.using.cutoff = remove.using.cutoff,
                                     percentage.dataset.removed = percentage.dataset.removed)

  colData(se)$Tumor_Molecular_Subtypes <- subtypes$consensusOV.subtypes
  # Print message to console
  cat("Classification complete: ", length(colData(se)$Tumor_Molecular_Subtypes ), "samples were classified for the consensus molecular subtypes.\n")

  return(se)
}
