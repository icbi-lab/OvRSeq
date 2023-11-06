#' Calculate Immunophenoscore (IPS) from a SummarizedExperiment
#'
#' This function calculates the Immunophenoscore (IPS) and its components scores
#' from gene expression data encapsulated in a `SummarizedExperiment` object.
#' It requires a specific set of genes and corresponding weights provided in a separate file.
#'
#' @param se A `SummarizedExperiment` object containing normalized gene expression data.
#' @return A data frame with samples as rows and calculated scores (WG, MHC, CP, EC, SC, MDSC, TREG, AZ, IPS) as columns.
#'
#' @importFrom SummarizedExperiment assay colData rowData
#' @references
#' Immunophenogram:
#' This R-script can be used to calculate Immunophenoscore (IPS) and generate Immunophenogram
#' from "EXPR.txt" and "IPS_genes.txt". The script and associated documentation can be found at
#' the following URL: https://github.com/icbi-lab/Immunophenogram
#' (C) ICBI, Medical University of Innsbruck, Biocenter, Division of Bioinformatics
#' Version 1.0 08.07.2016
#' @export
calculateIPS <- function(se) {

  data(IPS_genes)

  # Helper function to map raw scores to IPS
  ipsmap <- function(x) {
    if (x <= 0) return(0)
    if (x >= 3) return(10)
    round(x * 10 / 3, digits = 0)
  }

  # Read IPS genes and corresponding weights
  ips_genes <- IPS_genes

  # Prepare the expression matrix and gene list
  expression_data <- assay(se)
  sample_names <- colnames(expression_data)
  unique_ips_genes <- unique(ips_genes$GENE)

  # Initialize scores
  scores <- data.frame(matrix(ncol = 9, nrow = length(sample_names)))
  colnames(scores) <- c("MHC", "CP", "EC", "SC", "MDSC", "TREG", "AZ", "IPS")
  rownames(scores) <- sample_names

  # Calculate z-scores for each gene
  for (i in seq_along(sample_names)) {

    sample_expr <- expression_data[, i]
    z_scores <- (sample_expr - mean(sample_expr, na.rm = TRUE)) / sd(sample_expr, na.rm = TRUE)
    z_scores <- z_scores[ips_genes$GENE]
    z_scores_named <- setNames(z_scores, ips_genes$GENE)

    # Initialize vectors to store mean z-scores and weights
    MIG <- numeric(length(unique_ips_genes))
    WEIGHT <- numeric(length(unique_ips_genes))

    # Calculate mean z-scores and weights for each unique IPS gene
    for (j in seq_along(unique_ips_genes)) {
      gene_indices <- ips_genes$GENE[ips_genes$GENE == unique_ips_genes[j]]
      z_scores_subset <- z_scores_named[gene_indices]
      weights_subset <- ips_genes$WEIGHT[ips_genes$GENE == unique_ips_genes[j]]
      MIG[j] <- mean(z_scores_subset, na.rm = TRUE)
      WEIGHT[j] <- mean(weights_subset)
    }

    # Calculate weighted gene scores
    WG <- MIG * WEIGHT

    # Calculate component scores
    scores[i, "MHC"] <- mean(WG[1:10], na.rm = TRUE)
    scores[i, "CP"] <- mean(WG[11:20], na.rm = TRUE)
    scores[i, "EC"] <- mean(WG[21:24], na.rm = TRUE)
    scores[i, "SC"] <- mean(WG[25:26], na.rm = TRUE)
    scores[i, "MDSC"] <- WG[27] # Adjusted index for MDSC
    scores[i, "TREG"] <- WG[28] # Adjusted index for TREG

    # Calculate aggregate score
    scores[i, "AZ"] <- scores[i, "MHC"] + scores[i, "CP"] + scores[i, "EC"] + scores[i, "SC"]

    # Calculate final IPS score
    scores[i, "IPS"] <- ipsmap(scores[i, "AZ"])
  }

  # Save the results
  colData(se) <- cbind(colData(se), scores)

  # Check for NAs in IPS and calculate how many scores were computed
  num_calculated_ips <- sum(!is.na(scores$IPS))
  # Display a completion message with the count of non-NA IPS scores
  cat("IPS calculation complete: ", num_calculated_ips, "scores were calculated.\n")
  return(se)
}
