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

  # Wrapping each analysis block within a tryCatch to handle potential errors
  tryCatch({
    brcaness_classifier <- load_brcaness_classifier()
    brcaness_signature <- load_brcaness_signature()
    se <- classify_brcaness(se, brcaness_classifier, brcaness_signature)
  }, warning = function(w) {
    warning("Failed to classify BRCAness: ", w$message)
  })

  tryCatch({
    classifier_infiltration_status <- load_classifier_infiltration_status()
    immune_phenotype_signature <- load_tumor_immune_phenotype_signature()
    se <- classify_infiltration_status(se, classifier_infiltration_status, immune_phenotype_signature)
  }, warning = function(w) {
    warning("Failed to classify infiltration status: ", w$message)
  })

  tryCatch({
    se <- get_consensus_ov_subtypes(se, ids_type = "symbol")
  }, warning = function(w) {
    warning("Failed to determine consensus molecular subtypes: ", w$message)
  })

  tryCatch({
    se <- calculateIPS(se)
  }, warning = function(w) {
    warning("Failed to calculate IPS: ", w$message)
  })

  tryCatch({
    immune_signatures <- load_immune_signatures()
    se <- immune_signature_score(se, genelist = immune_signatures)
  }, warning = function(w) {
    warning("Failed to compute scores for standard immune gene signatures: ", w$message)
  })

  tryCatch({
    small_immune_signatures <- load_small_immune_signatures()
    se <- avg_expression_for_signature_se(se, small_immune_signatures)
  }, warning = function(w) {
    warning("Failed to compute average expression for small immune gene signatures: ", w$message)
  })

  # Perform immune deconvolution using multiple methods
  immunedeconv_methods <- c("quantiseq", "timer", "mcp_counter", "epic")
  for (method in immunedeconv_methods) {
    tryCatch({
      se <- deconvolute_immune(se, method)
    }, warning = function(w) {
      warning(paste("Failed immune deconvolution using", method, ":", w$message))
    })
  }

  # Specify the object names to be removed
  objs_to_remove <- c("brcaness_classifier", "classifier_infiltration",
                      "immune_signatures", "IPS_genes",
                      "brcaness_signature", "tumor_immune_phenotype_signature")

  # Remove objects if they exist in the global environment
  for (obj in objs_to_remove) {
    if (exists(obj, envir = .GlobalEnv)) {
      rm(list = obj, envir = .GlobalEnv)
    }
  }

  # Return the modified SummarizedExperiment object
  return(se)
}

