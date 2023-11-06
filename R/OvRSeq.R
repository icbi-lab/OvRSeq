

#se <- load_TCGA_OV()
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
  immunedeconv_methods <- c("quantiseq", "timer", "mcp_counter", "xcell", "epic")
  for (method in immunedeconv_methods) {
    se <- deconvolute_immune(se, method)
  }

  # Return the modified SummarizedExperiment object
  return(se)
}

