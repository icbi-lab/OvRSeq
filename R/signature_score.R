

#' Calculate enrichment scores for immune signatures
#'
#' This function calculates enrichment scores for a list of immune signatures
#' using the GSVA method. The user can provide their own gene list or use one
#' of the pre-defined ones. The function returns an updated SummarizedExperiment
#' object with the enrichment scores added to the colData.
#'
#' @param se A SummarizedExperiment object containing RNA sequencing data
#' @param method A character string indicating the method to use for calculating enrichment scores. Default is "ssgsea".
#' @param genelist A character vector with the gene list for the immune signature.
#'
#' @return A SummarizedExperiment object with the enrichment scores added to the colData.
#' @importFrom GSVA gsva
#' @export
immune_signature_score <- function(se, method = "ssgsea", genelist ){
  count_data <- assay(se)
  results <- gsva(count_data,
                 genelist,
                 method = 'ssgsea',
                 ssgsea.norm=F)
  # Create a success message
  success_message <- "Enrichment scores for immune signatures computed successfully."

  # Print the success message
  cat(success_message, "\n") # Using 'cat' to print the message followed by a newline character

  results <- t(results)
  colData(se) <- cbind(colData(se), results)
  return(se)
}

#' Compute average expression values for gene signatures and enrich colData of a SummarizedExperiment object
#'
#' Given a `SummarizedExperiment` object and a matrix or data frame of gene signatures,
#' this function computes the average expression values for each gene signature and enriches
#' the `colData` of the `SummarizedExperiment` object with the computed values.
#'
#' @param se A `SummarizedExperiment` object with gene expression data.
#' @param gene_sig A matrix or data frame where each column represents a gene signature and each cell contains a gene symbol.
#'        Empty cells should be represented as empty strings.
#'
#' @return A `SummarizedExperiment` object enriched with new columns in `colData`, each representing
#'         the average expression value of a given gene signature for each sample.
#'
#' @examples
#' # Assuming 'se' is a SummarizedExperiment object and 'gene_sig' is your gene signatures matrix
#' enriched_se <- avg_expression_for_signature_se(se, gene_sig, "MyGeneSet")
#'
#' @export
avg_expression_for_signature_se <- function(se, gmt) {

  # For each signature in the small_immune_signature list, compute average expression
  for (signature_name in names(gmt)) {
    gene_list <- gmt[[signature_name]]

    # Filter out empty entries if any exist
    gene_list <- gene_list[gene_list != ""]

    # Match rownames from se with gene_list to handle any missing genes
    matched_genes <- rownames(se)[rownames(se) %in% gene_list]

    # Get expression values for genes in the signature
    expression_genes <- assay(se)[matched_genes, , drop = FALSE]

    # Compute average expression across genes for each sample
    avg_expression <- colMeans(expression_genes, na.rm = TRUE)

    # Add this to colData of the SummarizedExperiment
    colData(se)[[paste0(signature_name,"_average")]] <- avg_expression
  }
  cat("The short immune signatures have been successfully computed and added to the dataset.\n")

  return(se)
}

#' Compute CYT to C1QA Ratio (C2C)
#'
#' This function calculates the cytolytic activity (CYT) to C1QA ratio (C2C)
#' for a given SummarizedExperiment object. The ratio is calculated using the
#' expression values of GZMB, PRF1, and C1QA. The C2C ratio is then added
#' to the `colData` of the provided SummarizedExperiment object.
#'
#' @param se A SummarizedExperiment object containing gene expression data.
#'           The object must contain expression data for GZMB, PRF1, and C1QA genes.
#'
#' @return Returns the SummarizedExperiment object with an additional column
#'         in `colData` named `ratio_CYT_C1QA`, representing the computed C2C ratio.
#'
#' @details The function computes the C2C ratio using the formula:
#'          C2C = 0.5 x (log2(TPM + 1 of GZMB) + log2(TPM + 1 of PRF1)) / log2(TPM + 1 of C1QA).
#'          It assumes that the expression data in the `SummarizedExperiment` object
#'          are either TPM (Transcripts Per Million) or intensity values suitable for
#'          the log transformation.
#'
#' @examples
#' # se is a SummarizedExperiment object with expression data for GZMB, PRF1, and C1QA
#' updatedSe <- computeCytC1qaRatio(se)
#'
#' @export
computeCytC1qaRatio <- function(se) {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment object.")
  }

  # Extract expression data for GZMB, PRF1, and C1QA
  exprData <- assay(se)
  requiredGenes <- c("GZMB", "PRF1", "C1QA")

  # Check if all required genes are present
  if (!all(requiredGenes %in% rownames(exprData))) {
    stop("Not all required genes (GZMB, PRF1, C1QA) are present in the assay data.")
  }

  # Calculate log2(TPM + 1) or log2 intensity values
  #logExprData <- log2(exprData[requiredGenes, ] + 1)

  # Compute the C2C ratio
  ratioCYT_C1QA <- 0.5 * (exprData["GZMB", ] + exprData["PRF1", ]) / exprData["C1QA", ]

  # Add the ratio to colData
  colData(se)$ratio_CYT_C1QA <- ratioCYT_C1QA

  cat("The CYT to C1QA Ratio have been successfully computed and added to the dataset.\n")

  return(se)
}

#' Compute Vulnerability Score for SummarizedExperiment Data
#'
#' This function calculates the vulnerability score for each sample in a `SummarizedExperiment` object. It requires the BRCAness probability and the ratio of cytolytic activity to C1QA (C2C ratio) as part of the dataset. If these are not present, the function computes them. The vulnerability score is computed using specified coefficients and a transformation function.
#'
#' @param se A `SummarizedExperiment` object that contains or can compute the BRCAness probability and the C2C ratio.
#'
#' @return The input `SummarizedExperiment` object with added columns for the mapped C2C ratio and the computed vulnerability score.
#'
#' @details
#' The function first checks if the `SummarizedExperiment` object contains `BRCAness_Prob` and `ratio_CYT_C1QA`. If not, it computes these using `classify_brcaness` and `computeCytC1qaRatio` functions, respectively. The vulnerability score is calculated using a logistic transformation of the C2C ratio and a linear combination of this transformed value with the BRCAness probability, using predefined coefficients. This score aims to reflect the patient's vulnerability based on their molecular profile.
#'
#' @examples
#' # se is a pre-loaded SummarizedExperiment object
#' se <- computeVulnerabilityScore(se)
#'
#' @importFrom SummarizedExperiment colData
#' @export
computeVulnerabilityScore <- function(se){
  # Check for required data and compute if necessary
  if (!("BRCAness_Prob" %in% colnames(colData(se)))) {
    #BRCAness probability is computed here and added to colData
    brcaness_classifier <- load_brcaness_classifier()
    brcaness_signature <- load_brcaness_signature()
    se <- classify_brcaness(se, brcaness_classifier, brcaness_signature)
  }
  if (!("ratio_CYT_C1QA" %in% colnames(colData(se)))) {
    #C2C ratio is computed here and added to colData
    se <- computeCytC1qaRatio(se)
  }

  # External functions
  mean_C2C <- 0.3013737
  sd_C2C <- 0.1359056

  mapping <- function(x) {
    v = x
    w = 1/(1+exp(-((x-mean_C2C)/(sd_C2C/pi))))
    return(w)
  }

  # Coefficients for vulnerability score calculation
  coefBRCAness <- 2.596728
  coefRatio <- 1.165859

  # Extract relevant data for plotting
  data <- colData(se)
  score <- (data$ratio_CYT_C1QA * coefRatio) + (data$BRCAness_Prob * coefBRCAness)
  colData(se)$mapped_ratio_CYT_C1QA <- mapping(data$ratio_CYT_C1QA)
  colData(se)$Vulnerability_Score <- score
  cat("The VulnerabilityScore have been successfully computed and added to the dataset.\n")

  return(se)
}
