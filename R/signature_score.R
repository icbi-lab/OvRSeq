

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

#' Compute Angiogenesis Score for SummarizedExperiment Data
#'
#' Calculates an angiogenesis score based on expression levels of seven angiogenic genes. This score ranges from 0 to 3 and is added to the `colData` of the `SummarizedExperiment` object.
#'
#' @param se `SummarizedExperiment` object containing the necessary gene expression data.
#' @param sample_id The identifier for the sample to compute the score.
#' @return `SummarizedExperiment` object with the new `angiogenesis_score` column added to its `colData`.
#'
#' @details
#' The function evaluates the expression of VEGFA, VEGFR2, PDGFA, PDGFB, PDGFRA, PDGFRB, and KIT.
#' It dichotomizes these based on a predefined threshold and calculates a score:
#' - 0 for 0 genes with high expression
#' - 1 for 1-3 genes with high expression
#' - 2 for 4-5 genes with high expression
#' - 3 for 6-7 genes with high expression
#' This score is indicative of the angiogenic potential of the tumor and has been linked to patient prognosis and response to angiogenesis inhibitors like bevacizumab.
#'
#' @references
#' Wieser V, Tsibulak I, Reimer DU, Zeimet AG, Fiegl H, Hackl H, Marth C. An angiogenic tumor phenotype
#' predicts poor prognosis in ovarian cancer. Gynecol Oncol. 2023 Mar;170:290-299.
#' doi: 10.1016/j.ygyno.2023.01.034. Epub 2023 Feb 7. PMID: 36758419.
#'
#' @examples
#' # se is a pre-loaded SummarizedExperiment object
#' sample_id <- "sample123"
#' se <- compute_angiogenesis_score(se, sample_id)
#'
#' @importFrom SummarizedExperiment colData assay
#' @export
compute_angiogenesis_score <- function(se) {
  required_genes <- c("VEGFA", "VEGFR2","KDR", "PDGFA", "PDGFB", "PDGFRA", "PDGFRB", "KIT")
  gene_thresholds <- list(
    VEGFA = 7.969714,
    VEGFR2 = 2.560207,
    KDR = 2.560207,
    PDGFA = 3.424363,
    PDGFB = 4.476959,
    PDGFRA = 3.695786,
    PDGFRB = 4.431066,
    KIT = 0.83184
  )

  # Extract the expression data for the required genes
  expr_data <- assay(se)
  required_genes <- intersect(rownames(expr_data), required_genes)
  if (length(required_genes) < 7 ) {
    warning("Some genes required for the angiogenesis signature are missing. The computed score may lack accuracy and should be interpreted with caution.")

  }
  expr_data <- expr_data[required_genes,]


  # Apply the score criteria based on expression levels
  number_of_high_exprs_gene <- apply(expr_data, 2, function(sample_expr) {
    sum(sapply(names(sample_expr), function(gene_name) {
      if (gene_name %in% names(gene_thresholds)) {
        sample_expr[gene_name] > gene_thresholds[[gene_name]]
      } else {
        warning(paste("Gene", gene_name, "not found in threshold list. Skipping this gene."))
        FALSE
      }
    }))
  })
  angiogenesis_score <- lapply(number_of_high_exprs_gene, function(x){
    if (x==0) {
      return(0)
    } else if (x > 0 && x < 4){
      return(1)
    } else if (x >= 4 && x < 6){
      return(2)
    } else if (x >= 6) {
      return(3)
    } else {
    return(NA)
    }
})


  # Add the angiogenesis score to colData
  colData(se)$angiogenesis_score <- unlist(angiogenesis_score)
  cat("The angiogenesis_score have been successfully computed and added to the dataset.\n")

  return(se)
}
