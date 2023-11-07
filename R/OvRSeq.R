#' Comprehensive RNA Sequencing-based Characterization of HGSOC Samples
#'
#' This function provides an all-encompassing analysis of high-grade serous ovarian cancer (HGSOC)
#' RNA-seq data. It integrates various predictive biomarkers and performs a multi-faceted immune
#' profiling, leveraging a 24-gene expression signature to classify BRCAness, subtype classification,
#' immune environment profiling, and immune cell deconvolution.
#'
#' @param se A SummarizedExperiment object containing RNA-seq data for HGSOC samples.
#' @param normalize A logical parameter indicating whether the assay data should be normalized.
#'                  If FALSE, the function will expect raw count data in integer form.
#' @return Returns a SummarizedExperiment object enriched with BRCAness classification,
#'         molecular subtype classification, immune phenotyping, IPS scoring, immune signatures,
#'         and immune cell deconvolution estimates.
#' @importFrom SummarizedExperiment assay colData
#' @examples
#' # Load data
#' se <- load_TCGA_OV()
#' # Run the OvRSeq analysis pipeline
#' se <- OvRSeq(se)
#'
#' @export
OvRSeq <- function(se, normalize = FALSE){
  # Validate input
  if (!inherits(se, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment object.")
  }

  # Warn about data normalization status
  assayData <- assay(se)
  dataType <- if (is.integer(assayData)) "raw counts" else "normalized"
  warning(sprintf("Data in 'se' appears to be %s.", dataType))

  if (is.numeric(assayData) && !normalize) {
    warning("Normalization is set to FALSE but data seems to be numeric. Ensure this is intentional.")
  }

  # Classify BRCAness status
  brcaness_classifier <- load_brcaness_classifier()
  brcaness_signature <- load_brcaness_signature()
  se <- classify_brcaness(se, brcaness_classifier, brcaness_signature)

  #Classify Infiltration status
  classifier_infiltration_status <- load_classifier_infiltration_status()
  immune_phenotype_signature <- load_tumor_immune_phenotype_signature()
  se <-classify_infiltration_status(se,classifier_infiltration_status,immune_phenotype_signature)

  # Determine consensus molecular subtypes
  se <- get_consensus_ov_subtypes(se, ids_type = "symbol")

  # Calculate the Immunophenoscore
  se <- calculateIPS(se)

  # Compute scores for standard immune gene signatures
  immune_signatures <- load_immune_signatures()
  se <- immune_signature_score(se, genelist = immune_signatures)

  # Compute average expression for small immune gene signatures
  small_immune_signatures <- load_small_immune_signatures()
  se <- avg_expression_for_signature_se(se, small_immune_signatures)

  # Perform immune deconvolution using multiple methods
  immunedeconv_methods <- c("quantiseq", "timer", "mcp_counter", "epic")
  for (method in immunedeconv_methods) {
    se <- deconvolute_immune(se, method)
  }

  # Return the modified SummarizedExperiment object
  return(se)
}

