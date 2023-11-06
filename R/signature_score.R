

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
                 ssgsea.norm=FALSE)
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
    colData(se)[[signature_name]] <- avg_expression
  }
  cat("The immune signatures have been successfully computed and added to the dataset.\n")

  return(se)
}


